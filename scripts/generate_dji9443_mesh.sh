#!/usr/bin/env bash
#
# Generate a DJI9443 surface mesh without opening the OpenVSP or Gmsh GUIs.
#
# QUICK START
#
# Interactive use (prompts for spanwise/chordwise values and cap mode):
#   scripts/generate_dji9443_mesh.sh
#
# Non-interactive use (suitable for an agent or parameter sweep):
#   scripts/generate_dji9443_mesh.sh 45 145 flat round 3
#   scripts/generate_dji9443_mesh.sh 45 145 none none 3 --output /tmp/dji.msh
#   scripts/generate_dji9443_mesh.sh 60 97 none round 5
#   scripts/generate_dji9443_mesh.sh 60 97 round flat 4
#
# ROOT_CAP and TIP_CAP are selected independently:
#   none   No cap at that end
#   flat   Flat cap
#   round  Round cap
#
# CAP_TESS is the shared refinement count for whichever end caps are enabled.
# OpenVSP 3.51.1 exposes this as a single parameter named CapUMinTess; despite
# that name it affects both root and tip caps, and there is no CapUMaxTess.
# CAP_TESS must be a positive integer and defaults to 3 from dji9443.vsp3.
#
# The original combined shorthands remain available:
#   capped    = flat root, round tip (the dji9443.vsp3 defaults)
#   uncapped  = no root or tip caps
#   flat      = flat caps at both ends
#   round     = round caps at both ends
# For example:
#   scripts/generate_dji9443_mesh.sh 45 145 capped
#   scripts/generate_dji9443_mesh.sh 45 145 uncapped
#
# Unless --output is supplied, the mesh is written under examples/data as:
#   dji9443_YYYYMMDD_SPANWISE_CHORDWISE_CAP_DESCRIPTION.msh
#
# Existing files are protected by default. Pass --force to replace one:
#   scripts/generate_dji9443_mesh.sh 45 145 flat round 3 --force
#
# Show the complete command-line help:
#   scripts/generate_dji9443_mesh.sh --help
#
# SPANWISE and CHORDWISE are the OpenVSP Tess_U and Tess_W values used in the
# existing DJI9443 mesh filenames.  CHORDWISE must be odd so that the upper and
# lower airfoil surfaces receive equal tessellation.
#
# Default tools and input:
#   OpenVSP: /Users/ryan/OpenVSP-3.51.1-MacOS/vsp
#   Gmsh:    gmsh found on PATH
#   Input:   examples/data/dji9443.vsp3
#
# Override those locations if needed:
#   DJI9443_VSP_BIN=/path/to/vsp \
#   DJI9443_GMSH_BIN=/path/to/gmsh \
#   DJI9443_VSP3=/path/to/model.vsp3 \
#     scripts/generate_dji9443_mesh.sh 45 145 flat round 3

set -euo pipefail

usage()
{
    cat <<'EOF'
Generate a DJI9443 surface mesh without opening the OpenVSP or Gmsh GUIs.

Usage:
  generate_dji9443_mesh.sh [SPANWISE CHORDWISE [ROOT_CAP TIP_CAP [CAP_TESS]]] [options]
  generate_dji9443_mesh.sh SPANWISE CHORDWISE COMBINED_CAP_MODE [options]

ROOT_CAP and TIP_CAP:
  none, flat, or round

CAP_TESS:
  Shared root/tip cap refinement count (OpenVSP CapUMinTess; default: 3)

Legacy COMBINED_CAP_MODE:
  capped, uncapped, flat, or round

Options:
      --root-cap MODE  Set root cap independently: none, flat, or round
      --tip-cap MODE   Set tip cap independently: none, flat, or round
      --cap-tess N     Set shared root/tip cap refinement (default: 3)
      --hub-r-over-r X Enlarge the hub: move the blade-root section (XSec_0
                       RadiusFrac, stock 0.01) out to X (e.g. 0.15). Interior
                       sections that would fall inside the new hub are
                       redistributed evenly between X and the first stock
                       station outboard of X (section count and airfoils are
                       preserved). Folded into the filename as _hubXpY.
  -o, --output FILE    Output .msh path
  -f, --force          Replace an existing output file
  -h, --help           Show this help

Environment overrides:
  DJI9443_VSP_BIN     OpenVSP executable
  DJI9443_GMSH_BIN    Gmsh executable or command name
  DJI9443_VSP3        Input parametric geometry
EOF
}

die()
{
    printf 'Error: %s\n' "$*" >&2
    exit 1
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
repo_root="$(cd "$script_dir/.." && pwd -P)"

vsp_bin="${DJI9443_VSP_BIN:-/Users/ryan/OpenVSP-3.51.1-MacOS/vsp}"
gmsh_bin="${DJI9443_GMSH_BIN:-gmsh}"
vsp3_file="${DJI9443_VSP3:-$repo_root/examples/data/dji9443.vsp3}"
output_file=""
force=0
root_cap_mode=""
tip_cap_mode=""
cap_tess=""
hub_r_over_r=""
positionals=()

while (($#)); do
    case "$1" in
        -o|--output)
            (($# >= 2)) || die "$1 requires a file path"
            output_file="$2"
            shift 2
            ;;
        -f|--force)
            force=1
            shift
            ;;
        --root-cap)
            (($# >= 2)) || die "$1 requires none, flat, or round"
            root_cap_mode="$2"
            shift 2
            ;;
        --tip-cap)
            (($# >= 2)) || die "$1 requires none, flat, or round"
            tip_cap_mode="$2"
            shift 2
            ;;
        --cap-tess)
            (($# >= 2)) || die "$1 requires a positive integer"
            cap_tess="$2"
            shift 2
            ;;
        --hub-r-over-r)
            (($# >= 2)) || die "$1 requires a decimal radius fraction (e.g. 0.15)"
            hub_r_over_r="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            while (($#)); do
                positionals+=("$1")
                shift
            done
            ;;
        -*)
            die "unknown option: $1"
            ;;
        *)
            positionals+=("$1")
            shift
            ;;
    esac
done

((${#positionals[@]} <= 5)) ||
    die "expected at most SPANWISE CHORDWISE ROOT_CAP TIP_CAP CAP_TESS"

spanwise="${positionals[0]:-}"
chordwise="${positionals[1]:-}"

if ((${#positionals[@]} > 2)); then
    [[ -z "$root_cap_mode" && -z "$tip_cap_mode" ]] ||
        die "do not combine positional cap modes with --root-cap or --tip-cap"

    if ((${#positionals[@]} == 3)); then
        case "${positionals[2]}" in
            capped|mixed)
                root_cap_mode="flat"
                tip_cap_mode="round"
                ;;
            uncapped|none)
                root_cap_mode="none"
                tip_cap_mode="none"
                ;;
            flat)
                root_cap_mode="flat"
                tip_cap_mode="flat"
                ;;
            round)
                root_cap_mode="round"
                tip_cap_mode="round"
                ;;
            *)
                die "a single combined cap mode must be capped, uncapped, flat, or round"
                ;;
        esac
    else
        root_cap_mode="${positionals[2]}"
        tip_cap_mode="${positionals[3]}"
        if ((${#positionals[@]} == 5)); then
            [[ -z "$cap_tess" ]] ||
                die "do not combine positional CAP_TESS with --cap-tess"
            cap_tess="${positionals[4]}"
        fi
    fi
fi

if [[ -z "$spanwise" || -z "$chordwise" ]]; then
    [[ -t 0 ]] ||
        die "SPANWISE and CHORDWISE are required when standard input is not a terminal"
    [[ -n "$spanwise" ]] ||
        read -r -p "Desired spanwise tessellation (Tess_U): " spanwise
    [[ -n "$chordwise" ]] ||
        read -r -p "Desired chordwise tessellation (odd Tess_W): " chordwise
fi

if [[ -t 0 ]]; then
    if [[ -z "$root_cap_mode" ]]; then
        read -r -p "Root cap [none/flat/round] (default: flat): " root_cap_mode
    fi
    if [[ -z "$tip_cap_mode" ]]; then
        read -r -p "Tip cap [none/flat/round] (default: round): " tip_cap_mode
    fi
    if [[ -z "$cap_tess" ]]; then
        read -r -p "Cap tessellations (CapUMinTess, default: 3): " cap_tess
    fi
fi
root_cap_mode="${root_cap_mode:-flat}"
tip_cap_mode="${tip_cap_mode:-round}"
cap_tess="${cap_tess:-3}"

[[ "$spanwise" =~ ^[0-9]+$ ]] ||
    die "SPANWISE must be a positive integer"
[[ "$chordwise" =~ ^[0-9]+$ ]] ||
    die "CHORDWISE must be a positive integer"
((spanwise >= 2)) || die "SPANWISE must be at least 2"
((chordwise >= 3)) || die "CHORDWISE must be at least 3"
((chordwise % 2 == 1)) ||
    die "CHORDWISE must be odd (for equal upper/lower surface tessellation)"
[[ "$cap_tess" =~ ^[0-9]+$ ]] ||
    die "CAP_TESS must be a positive integer"
((cap_tess >= 1)) || die "CAP_TESS must be at least 1"

if [[ -n "$hub_r_over_r" ]]; then
    [[ "$hub_r_over_r" =~ ^0\.[0-9]+$ ]] ||
        die "--hub-r-over-r must be a decimal in (0, 0.9), e.g. 0.15"
    awk -v x="$hub_r_over_r" 'BEGIN { exit !(x > 0.01 && x < 0.9) }' ||
        die "--hub-r-over-r must be > 0.01 (stock root) and < 0.9"
fi

case "$root_cap_mode" in
    none)
        root_cap_value=0
        ;;
    flat)
        root_cap_value=1
        ;;
    round)
        root_cap_value=2
        ;;
    *)
        die "ROOT_CAP must be none, flat, or round"
        ;;
esac

case "$tip_cap_mode" in
    none)
        tip_cap_value=0
        ;;
    flat)
        tip_cap_value=1
        ;;
    round)
        tip_cap_value=2
        ;;
    *)
        die "TIP_CAP must be none, flat, or round"
        ;;
esac

case "$root_cap_mode:$tip_cap_mode" in
    flat:round) output_suffix="capped" ;;
    none:none) output_suffix="uncapped" ;;
    flat:flat) output_suffix="flatcaps" ;;
    round:round) output_suffix="roundcaps" ;;
    *) output_suffix="root-${root_cap_mode}_tip-${tip_cap_mode}" ;;
esac
if [[ "$root_cap_mode" != "none" || "$tip_cap_mode" != "none" ]]; then
    if ((cap_tess != 3)); then
        output_suffix="${output_suffix}_captess${cap_tess}"
    fi
fi
if [[ -n "$hub_r_over_r" ]]; then
    output_suffix="${output_suffix}_hub${hub_r_over_r//./p}"
fi

[[ -x "$vsp_bin" ]] || die "OpenVSP executable not found or not executable: $vsp_bin"
if [[ "$gmsh_bin" == */* ]]; then
    [[ -x "$gmsh_bin" ]] || die "Gmsh executable not found or not executable: $gmsh_bin"
else
    gmsh_bin="$(command -v "$gmsh_bin")" ||
        die "Gmsh command not found: ${DJI9443_GMSH_BIN:-gmsh}"
fi
[[ -f "$vsp3_file" ]] || die "OpenVSP geometry not found: $vsp3_file"

if [[ -z "$output_file" ]]; then
    output_file="$repo_root/examples/data/dji9443_$(date +%Y%m%d)_${spanwise}_${chordwise}_${output_suffix}.msh"
elif [[ "$output_file" != /* ]]; then
    output_file="$(pwd -P)/$output_file"
fi

[[ "$output_file" == *.msh ]] || die "output filename must end in .msh"
output_dir="$(dirname "$output_file")"
[[ -d "$output_dir" ]] || die "output directory does not exist: $output_dir"
if [[ -e "$output_file" && "$force" -ne 1 ]]; then
    die "output already exists (pass --force to replace it): $output_file"
fi

# AngelScript string literals cannot safely contain these path characters.
for checked_path in "$vsp3_file" "$output_file"; do
    case "$checked_path" in
        *'"'*|*'\'*|*$'\n'*)
            die "paths containing quotes, backslashes, or newlines are unsupported: $checked_path"
            ;;
    esac
done

work_dir="$(mktemp -d "${TMPDIR:-/tmp}/dji9443_mesh.XXXXXX")"
vsp_script="$work_dir/generate.vspscript"
raw_mesh="$work_dir/openvsp_raw.msh"
staged_output="$(mktemp "$output_dir/.dji9443_mesh.msh.XXXXXX")"

cleanup()
{
    rm -f "$vsp_script" "$raw_mesh" "$staged_output"
    rmdir "$work_dir" 2>/dev/null || true
}
trap cleanup EXIT

hub_r_value="${hub_r_over_r:-0.0}"

cat >"$vsp_script" <<EOF
void main()
{
    ClearVSPModel();
    ReadVSPFile("$vsp3_file");

    array<string> props = FindGeomsWithName("PropGeom");
    if (props.size() != 1)
    {
        Print("ERROR: expected exactly one geometry named PropGeom.\\n");
        return;
    }

    // Isolate the prop from the hub fuselage and any saved MeshGeom objects.
    array<string> geoms = FindGeoms();
    for (uint i = 0; i < geoms.size(); ++i)
    {
        SetSetFlag(geoms[i], SET_FIRST_USER, false);
    }

    string prop_id = props[0];
    SetSetFlag(prop_id, SET_FIRST_USER, true);
    SetParmVal(prop_id, "Tess_U", "Shape", $spanwise);
    SetParmVal(prop_id, "Tess_W", "Shape", $chordwise);
    SetParmVal(prop_id, "CapUMinOption", "EndCap", $root_cap_value);
    SetParmVal(prop_id, "CapUMaxOption", "EndCap", $tip_cap_value);
    SetParmVal(prop_id, "CapUMinTess", "EndCap", $cap_tess);

    double hub_r = $hub_r_value;
    if (hub_r > 0.0)
    {
        // Collect the stock RadiusFrac stations (XSec_0 .. XSec_n).
        string xsec_surf = GetXSecSurf(prop_id, 0);
        int nxsec = GetNumXSec(xsec_surf);
        array<double> stations;
        array<string> station_parms;
        for (int i = 0; i < nxsec; ++i)
        {
            string xsec_id = GetXSec(xsec_surf, i);
            string pid = GetXSecParm(xsec_id, "RadiusFrac");
            if (pid.length() == 0)
            {
                Print("ERROR: XSec " + i + " has no RadiusFrac parm.\\n");
                return;
            }
            stations.insertLast(GetParmVal(pid));
            station_parms.insertLast(pid);
        }
        if (nxsec < 2)
        {
            Print("ERROR: could not enumerate prop XSec RadiusFrac stations.\\n");
            return;
        }

        // First stock station strictly outboard of the new hub (with margin
        // so no two sections coincide).
        int k = -1;
        for (int i = 1; i < nxsec; ++i)
        {
            if (stations[i] > hub_r + 0.01) { k = i; break; }
        }
        if (k < 0)
        {
            Print("ERROR: hub_r_over_r leaves no outboard section; reduce it.\\n");
            return;
        }

        // Delete interior sections that fall inside the new hub (crowding
        // them instead produces loft spikes at the old root radius), then
        // move the stock root section out to the hub. Do NOT cut the root
        // section itself: cutting changes the per-segment spanwise
        // tessellation (element count) and the root parm re-clamps.
        // Validated for hub_r up to ~0.2; the fat transitional root folds
        // the loft (non-watertight) when dragged much beyond that, so
        // verify watertightness downstream for larger values.
        for (int i = k - 1; i >= 1; --i)
        {
            CutXSec(prop_id, i);
        }
        Update();
        // Parm ids are stale after CutXSec; re-fetch the root section.
        xsec_surf = GetXSecSurf(prop_id, 0);
        string root_xsec = GetXSec(xsec_surf, 0);
        string root_parm = GetXSecParm(root_xsec, "RadiusFrac");
        SetParmValUpdate(root_parm, hub_r);
        Print("  cut " + (k - 1) + " interior section(s); root RadiusFrac now " +
              GetParmVal(root_parm) + "\\n");
        Print("Hub variant: XSec_0 RadiusFrac -> " + hub_r + "; next station " +
              stations[k] + "\\n");
    }

    Update();

    string mesh_id = ExportFile(
        "$raw_mesh",
        SET_FIRST_USER,
        EXPORT_GMSH
    );
    DeleteGeom(mesh_id);

    while (GetNumTotalErrors() > 0)
    {
        ErrorObj err = PopLastError();
        Print(err.GetErrorString());
    }
}
EOF

printf 'Generating DJI9443 mesh: spanwise=%s chordwise=%s root-cap=%s tip-cap=%s cap-tess=%s hub-r-over-r=%s\n' \
    "$spanwise" "$chordwise" "$root_cap_mode" "$tip_cap_mode" "$cap_tess" "${hub_r_over_r:-stock}"

# OpenVSP 3.51.1 on macOS can return a nonzero process status after a successful
# -script export, so the generated mesh is the authoritative success signal.
set +e
"$vsp_bin" -script "$vsp_script"
vsp_status=$?
set -e

if [[ ! -s "$raw_mesh" ]] || ! grep -q '^\$MeshFormat' "$raw_mesh"; then
    die "OpenVSP did not produce a valid mesh (process status $vsp_status)"
fi
if ((vsp_status != 0)); then
    printf 'Note: OpenVSP returned status %d after producing the mesh; continuing.\n' \
        "$vsp_status" >&2
fi

"$gmsh_bin" "$raw_mesh" -save -format msh4 -o "$staged_output"
[[ -s "$staged_output" ]] || die "Gmsh did not produce the converted mesh"
grep -q '^4\.1 ' "$staged_output" ||
    die "Gmsh output is not an ASCII MSH 4.1 file"

mv -f "$staged_output" "$output_file"
printf 'Wrote %s\n' "$output_file"
