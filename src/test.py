from google import genai
from google.genai.types import HttpOptions

# Explicitly routing to the beta endpoint where 3.1 lives
client = genai.Client(
    vertexai=True,
    project="305675469459",
    location="us-central1",
    http_options=HttpOptions(api_version="v1beta1")
)

print("Connecting to Gemini 3.1 Pro Preview...")
response = client.models.generate_content(
    model="gemini-3.1-pro-preview",
    contents="Hello! Acknowledge connection."
)

print(f"Success: {response.text}")
