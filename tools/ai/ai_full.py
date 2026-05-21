import os
import time
from google import genai
from google.genai import types

# --- 2026 PRICING REGISTRY (Implicit Cache Discount Rates) ---
PRICING = {
    "gemini-3.1-pro-preview": {"in": 0.40, "out": 18.00},
    "gemini-3-flash-preview": {"in": 0.20, "out": 3.00},
    "gemini-2.5-flash": {"in": 0.10, "out": 0.40},
    "gemini-3.1-flash-lite-preview": {"in": 0.01, "out": 0.10}
}

PROJECT_ID = "305675469459"
LOCATION = "global"
#LOCATION = "us-central1"
client = genai.Client(
    vertexai=True,
    project=PROJECT_ID,
    location=LOCATION,
)

def load_code(directory="./"):
    """Crawls the directory to build a single codebase string."""
    code_parts = []
    exts = ('.cpp', '.h', '.cu', '.py', '.jl', '.f90')
    for root, _, files in os.walk(directory):
        if any(x in root for x in ['build', '.git', '__pycache__']): continue
        for f in files:
            if f.endswith(exts):
                path = os.path.relpath(os.path.join(root, f), directory)
                try:
                    with open(os.path.join(root, f), 'r', encoding='utf-8') as s:
                        code_parts.append(f"--- FILE: {path} ---\n{s.read()}\n")
                except Exception as e:
                    print(f"Skipping {path}: {e}")
    return "\n".join(code_parts)

def create_chat(model_id, thinking="MEDIUM", codebase=""):
    """Spins up a chat session with the codebase embedded in the system prompt."""
    print(f"🔄 Booting {model_id}...")

    # Base configuration with your PhD Persona
    config_kwargs = {
        "system_instruction": f"You are a PhD expert in computational aerodynamics. Here is my codebase:\n{codebase}",
        "temperature": 0.2  # Low temperature for precise code/math
    }

    # Only apply Thinking Levels to the Gemini 3/3.1 families
    if "3.1" in model_id or "3-" in model_id:
        config_kwargs["thinking_config"] = types.ThinkingConfig(thinking_level=thinking)

    config = types.GenerateContentConfig(**config_kwargs)
    return client.chats.create(model=model_id, config=config), model_id

# --- MASTER LOOP ---
print("🚀 Master CFD Terminal (Plan B: Implicit Caching)")
print("Options: 'model [pro|flash3|flash2|lite]' | 'think [low|medium|high]' | 'clear' | 'exit'")

print("📦 Reading codebase from local disk...")
global_codebase = load_code()

# Start with 3 Flash Preview as the speedy default
chat, current_model = create_chat("gemini-2.5-flash", "MEDIUM", global_codebase)
current_think = "MEDIUM"

while True:
    cmd = input(f"\n({current_model} | {current_think}) > ").strip()

    # 1. Handle Model Switching
    if cmd.lower().startswith("model"):
        target = cmd.split()[-1]
        m_map = {
            "pro": "gemini-3.1-pro-preview",
            "flash3": "gemini-3-flash-preview",
            "flash2": "gemini-2.5-flash",
            "lite": "gemini-3.1-flash-lite-preview"
        }
        if target in m_map:
            chat, current_model = create_chat(m_map[target], current_think, global_codebase)
        else:
            print(f"⚠️ Unknown model. Available options: {list(m_map.keys())}")
        continue

    # 2. Handle Thinking Level Adjustments
    if cmd.lower().startswith("think"):
        target_think = cmd.split()[-1].upper()
        if target_think not in ["LOW", "MEDIUM", "HIGH"]:
            print("⚠️ Invalid thinking level. Choose LOW, MEDIUM, or HIGH.")
            continue

        if "3.1" not in current_model and "3-" not in current_model:
            print("⚠️ Thinking levels are only supported on Gemini 3 models. Switch to 'pro', 'flash3', or 'lite' first.")
            continue

        current_think = target_think
        chat, current_model = create_chat(current_model, current_think, global_codebase)
        print(f"🧠 Reasoning depth set to {current_think}")
        continue

    # 3. Handle Chat Management
    if cmd.lower() == 'clear':
        print("🗑️ Wiping chat history. Codebase remains loaded via System Prompt.")
        chat, current_model = create_chat(current_model, current_think, global_codebase)
        continue

    if cmd.lower() in ['exit', 'quit']:
        print("👋 Session ended. Zero rent incurred today.")
        break

    if not cmd:
        continue

    # 4. Send the message and calculate costs
    try:
        response = chat.send_message(cmd)

        # Pull pricing based on the active model
        p = PRICING.get(current_model, {"in": 0.0, "out": 0.0})

        # Calculate tokens
        prompt_tokens = response.usage_metadata.prompt_token_count
        out_tokens = response.usage_metadata.candidates_token_count

        in_cost = (prompt_tokens / 1_000_000) * p["in"]
        out_cost = (out_tokens / 1_000_000) * p["out"]

        print(f"\nGemini: {response.text}")
        print(f"\n[💰 Turn Cost: ${in_cost+out_cost:.4f} | Tokens: {prompt_tokens:,} in / {out_tokens:,} out]")

    except Exception as e:
        print(f"\n❌ API Error: {e}")
