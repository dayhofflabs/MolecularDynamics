# Github Actions Workflow

The files in this folder control the tests that we run in GHA.  We'd like to use the same base image for our tests as we use locally for the dev container and on Batch for any job containers, but we can't because that image is too big to fit in the GHA runners. As a workaround, we build a minimal image with whatever packages are needed to make it work. Those packages are pip-installed, and represent the bare minimum to make the tests pass. To test that locally, run the script `.devcontainer/scripts/3_test_pytest_in_github_actions_container.sh`.
