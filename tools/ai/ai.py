import vertexai
from vertexai.generative_models import GenerativeModel

# Initialize with your project info
vertexai.init(project="project-1-489216", location="us-central1")

# Use Gemini 3 Pro for complex aerodynamic logic
model = GenerativeModel("gemini-2.5-flash")
# model = GenerativeModel("gemini-3.1-pro-preview")

prompt = """
Tell my about your day.
"""

response = model.generate_content(prompt)
print(response.text)
