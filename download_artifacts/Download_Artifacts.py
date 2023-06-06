import os
import io
import zipfile
import requests
from github import Github

def download_artifacts(token, repo_owner, repo_name, workflow_run_id):
    # Initialize GitHub API client
    gh = Github(token)
    repo = gh.get_repo(f"{repo_owner}/{repo_name}")

    # Get the specified workflow run artifacts
    workflow_run = repo.get_workflow_run(workflow_run_id)
    artifacts = workflow_run.get_artifacts()

    # Get the current file path and set it to the Data_downloaded folder
    file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "Data_downloaded")
    os.makedirs(file_path, exist_ok=True)

    # Download and extract all artifacts
    for artifact in artifacts:
        response = requests.get(artifact.archive_download_url, headers={"Authorization": f"token {token}"})
        
        if response.status_code == 200:
            with zipfile.ZipFile(io.BytesIO(response.content)) as zf:
                zf.extractall(file_path)
                print(f"Downloaded and extracted {artifact.name}")
        else:
            print(f"Error downloading {artifact.name}: {response.status_code}")

if __name__ == "__main__":
    GITHUB_TOKEN = "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"
    REPO_OWNER = "OptHuang"
    REPO_NAME = "option-pricing-regime-jump_matlab"
    WORKFLOW_RUN_ID = 5148212386

    download_artifacts(GITHUB_TOKEN, REPO_OWNER, REPO_NAME, WORKFLOW_RUN_ID)
