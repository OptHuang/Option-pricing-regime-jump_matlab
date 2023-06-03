import os
import requests
from github import Github

# Set your GitHub token and repository details
GITHUB_TOKEN = "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"
REPO_OWNER = "OptHuang"
REPO_NAME = "Option-pricing-regime-jump_matlab"
WORKFLOW_RUN_ID = 5148212386  # Replace with the desired workflow run ID

# Initialize GitHub API client
gh = Github(GITHUB_TOKEN)
repo = gh.get_repo(f"{REPO_OWNER}/{REPO_NAME}")

# Get the specified workflow run artifacts
workflow_run = repo.get_workflow_run(WORKFLOW_RUN_ID)
artifacts = workflow_run.get_artifacts()

# Get the current file path
file_path = os.path.dirname(os.path.abspath(__file__))

# Download all artifacts
for artifact in artifacts:
    response = requests.get(artifact.archive_download_url, headers={"Authorization": f"token {GITHUB_TOKEN}"})
    
    # Save the .mat file as binary
    with open(os.path.join(file_path, artifact.name), "wb") as f:
        f.write(response.content)

    print(f"Downloaded {artifact.name}")
