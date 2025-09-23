# Contributing to Goombay

Thank you for your interest in contributing to Goombay! We welcome contributions,
whether you're fixing a bug, improving documentation, or adding new features.
This guide will help you get started with contributing to the project.

We recommend using [uv](https://docs.astral.sh/uv/getting-started/installation/)
for managing dependencies and running the project locally.

After installing uv, check out [uv commands](https://docs.astral.sh/uv/reference/cli/).
The most important will be `uv run` to run project files, `uv sync` to sync dependencies,
`uv add` to add dependencies, `uv remove` to remove dependencies, and `uv pip`.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [How to Contribute](#how-to-contribute)
- [Bug Reports](#bug-reports)
- [Feature Requests](#feature-requests)
- [Pull Requests](#pull-requests)
- [Running Tests](#running-tests)
- [Code Style](#code-style)
- [License](#license)
- [Community](#community)

## Code of Conduct

By participating in this project, you agree to abide by our
[Code of Conduct](https://github.com/lignum-vitae/goombay/blob/master/docs/CODE_OF_CONDUCT.md).
Please take a moment to familiarize yourself with it.

## How to Contribute

### Bug Reports

If you find a bug or unexpected behavior, please open an issue in the GitHub
repository. When reporting a bug, provide the following information:

- A description of the problem.
- Steps to reproduce the issue (if applicable).
- Any relevant error messages.
- The version of the library you're using.
- The Python version you're using.

### Feature Requests

If you have an idea for a new feature or enhancement, please open an issue
describing the feature and why you think it would be useful.
We encourage open discussions before starting to code a new feature.

### Pull Requests

To contribute code:

#### 1. Open a new Issue following the above-mentioned guidelines

#### 2. Fork the repository to your own GitHub account

#### 3. Clone your fork locally

```nginx
git clone https://github.com/YOUR_USERNAME/goombay.git
```

#### 4. Keep your fork up to date with the main repository

```nginx
# Add the main Goombay repository as the upstream branch (to pull future any updates from)
git remote add upstream https://github.com/lignum-vitae/goombay.git
# Get latest changes
git fetch upstream
# Verify your remotes
git remote
```

#### 5. Checkout your branch for changes

##### Choose ONE of the following commands

```nginx
# Creates a new branch that stays in sync with the develop branch of the main repository
git checkout -b <feature-name> upstream/develop

# Checks out existing branch if you already have a branch locally
git checkout <feature-name>
```

#### 6. Make your changes in your local repository and run in your local environment

```nginx
# Downloads project as editable, which allows local imports. Run this command from the root directory.
uv pip install -e ".[dev]"

# Ensure dependencies match the project's using uv sync
uv sync --all-extras

# Files may be run preferentially run using uv
uv run editdistance.py

# Or as a module
python -m goombay.algorithms.editdistance

# Or as a script
python editdistance.py
```

#### 7. Test your changes. Make sure all tests pass. (See [Running Tests](#running-tests))

#### 8. Commit your changes with a descriptive commit message

```nginx
# Gets latest changes from main biobase project if you've set up an upstream branch as detailed above
git fetch upstream
# We recommend individually adding each file with modifications
git add <filename>
# Commit files after all files with modifications have been added
git commit -m "Add feature: description of change"`
```

#### 🚨 Using git add .

While `git add .` is convenient for adding all modified files, it can lead
to messy commits. Consider using it only when:

- You've reviewed all changes
- You're certain about each modification
- You've checked git status first
- Your .gitignore is properly configured

#### 9. Rebase Your Development Branch with the Latest Upstream changes

```nginx
# Make sure all is committed (or stashed) as necessary on this branch
git rebase -i upstream/develop feature-name
```

You may need to resolve conflicts that occur when both a file on the development
trunk and one of the files in your branch have been changed.
Edit each conflicting file to resolve the differences, then continue the rebase.
Each file will need to be "added" to mark the conflict as resolved:

```nginx
# Resolve conflicts in each file, then:
git add <resolved-filename>
git rebase --continue
```

The behavior of `git pull` can be adjusted to deal with merge conflicts by
rebasing them instead of merging them by running the command
`git config pull.rebase true` from the base of this repo.

This will only change the behavior for this repo. To change it for all your current
repositories, run the command `git config --global pull.rebase true`.

This change can be confirmed by checking that this is in your `.git/config` file:

```toml
[pull]
	rebase = true
```

#### 10. Push your branch to your fork on GitHub

```nginx
git push -f origin feature-name
```

#### 11. Open a Pull Request (PR) from your branch to the develop branch of the original Goombay repository.

- You may need to click the `compare across forks` link under the
  `Compare changes` header that populates
  when you click `New pull request` to see your local repo fork.

#### 12. In your PR description, include

- A summary of the changes.
- Any relevant issue numbers (e.g., fixes #123).
- Information about tests and validation.

## Running Tests

To ensure your changes work correctly, you can run the tests before
submitting a PR.

Install development dependencies:

```nginx
uv sync --all-extras
```
Run all tests using one of the following commands from the root directory:

```nginx
python -m unittest discover tests
```

or if you're using py:

```nginx
py -m unittest discover tests
```

or if using pytest:

```nginx
pytest tests\
```

This will run the tests in the tests directory.

## Code Style

We use [**Black**](https://black.readthedocs.io/en/stable/) to automatically
format Python code. Black enforces a consistent style that may differ slightly
from PEP 8, and its formatting decisions are not configurable.

Please run Black on your code before submitting a pull request. You can install
it via pip:

```shell
# If black is not already installed, pip install black
pip install black

# Run this command from the root in order to format all python files
black .
```

Additional guidelines:

- Use meaningful variable and function names.
- Ensure all new features are covered by tests.
- Update documentation if your changes affect usage or the API.

## License

By contributing to Goombay, you agree that your contributions will be licensed
under the MIT License, as outlined in the
[LICENSE](https://github.com/lignum-vitae/goombay/blob/master/LICENSE) file.

## Community

We encourage contributions from everyone, and we strive to maintain a welcoming
and inclusive community. If you have any questions, need help,
or want to discuss ideas, feel free to reach out via issues or the repository discussions.

Thank you for contributing to Goombay! Your help improves the project for everyone!
