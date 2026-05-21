import os
from google import genai

# --- CONFIGURATION ---
# Replace with your actual Project ID from the Vertex AI console
PROJECT_ID = "305675469459"
LOCATION = "global" # Or your preferred region
MODEL_ID = "gemini-3.1-pro-preview"

def run_terminal_chat():
    # Initialize the client for Vertex AI
    client = genai.Client(
        vertexai=True,
        project=PROJECT_ID,
        location=LOCATION
    )

    # Start a stateful chat session
    chat_session = client.chats.create(model=MODEL_ID)

    print(f"\n[ Connected to {MODEL_ID} ]")
    print("Type 'exit' or 'quit' to end the session.\n")

    while True:
        try:
            user_input = input("You: ").strip()

            if not user_input:
                continue
            if user_input.lower() in ["exit", "quit"]:
                print("Ending chat. Goodbye!")
                break

            print("Gemini: ", end="", flush=True)

            # Use the dedicated stream method instead of a keyword argument
            for chunk in chat_session.send_message_stream(user_input):
                if chunk.text:
                    print(chunk.text, end="", flush=True)

            print("\n")
        except KeyboardInterrupt:
            print("\nSession interrupted.")
            break
        except Exception as e:
            print(f"\nAn error occurred: {e}")

if __name__ == "__main__":
    run_terminal_chat()
