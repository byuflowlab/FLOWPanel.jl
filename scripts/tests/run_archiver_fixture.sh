#!/usr/bin/env bash
set -e
rm -rf fx && mkdir -p fx/data fx/arch fx/BR && cd fx
printf 'name = "FLOWPanel"\nuuid = "6be8c882-484d-4309-b349-d23112750151"\n' > Project.toml
mkrun() {
  r=$1; n=$2; d=data/$r
  mkdir -p "$d/${r}_body1" "$d/${r}_wake1" "$d/${r}_wake1_particles" "$d/monitors"
  i=1; while [ "$i" -le "$n" ]; do
    printf 'body %s step %s padding\n' "$r" "$i" > "$d/${r}_body1/${r}_body1.$i.vtu"
    printf 'wake1 %s %s\n' "$r" "$i" > "$d/${r}_wake1/${r}_wake1.1.$i.vts"
    printf 'wake2 %s %s\n' "$r" "$i" > "$d/${r}_wake1/${r}_wake1.2.$i.vts"
    printf 'part %s %s\n' "$r" "$i" > "$d/${r}_wake1_particles/${r}_wake1_particles.$i.vtp"
    printf 'idx %s %s\n'  "$r" "$i" > "$d/${r}_wake1.$i.vtm"
    i=$((i+1))
  done
  printf 'pvd\n' > "$d/${r}_body1.pvd"
  printf 'rev,CT\n1,0.05\n' > "$d/${r}_CT_vs_rev.csv"
  printf 'x=1\n' > "$d/${r}.metadata.toml"
  printf 'x=1\n' > "$d/${r}_case_metadata.toml"
  printf 'mon\n' > "$d/monitors/${r}_force.csv"
  printf 'monvtk\n' > "$d/monitors/${r}_never_touch.vtu"
}
mkrun finished_a 8; mkrun finished_b 3; mkrun protected_c 3; mkrun live_d 3; mkrun recent_f 7; mkrun recent_g 7
mkdir -p data/novtk_e/monitors
printf 'rev,CT\n1,0.05\n' > data/novtk_e/novtk_e_CT_vs_rev.csv
find data -exec touch -t 202601010000 {} + 2>/dev/null
# Ages chosen to straddle the two thresholds: <2 h is RECENT-HOT (never
# approvable), 2-24 h is RECENT (approvable by name), >=24 h archives itself.
ago() { date -v-"$1"H +%Y%m%d%H%M 2>/dev/null || date -d "$1 hours ago" +%Y%m%d%H%M; }
touch data/live_d/live_d_body1/live_d_body1.3.vtu                     # now  -> RECENT-HOT
touch -t "$(ago 5)" data/recent_f/recent_f_body1/recent_f_body1.7.vtu  # 5 h  -> RECENT
touch -t "$(ago 6)" data/recent_g/recent_g_body1/recent_g_body1.7.vtu  # 6 h  -> RECENT
printf '# protected\nprotected_c\n' > BR/protect.txt
cp -R data ref

# A SECOND checkout, holding a run with the SAME NAME as one in the first --
# the collision the per-checkout slug directories have to keep apart.
mkdir -p ../fx2/data && cd ../fx2 && printf 'name = "FLOWPanel"\nuuid = "6be8c882-484d-4309-b349-d23112750151"\n' > Project.toml
mkrun finished_a 6
mkrun other_only 7
find data -exec touch -t 202601010000 {} + 2>/dev/null
cd ../fx
