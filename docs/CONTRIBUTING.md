# Contributing to Goombay
Thank you for your interest in contributing to Goombay! We welcome contributions, whether you're fixing a bug, improving documentation, or adding new features. This guide will help you get started with contributing to the project.

## Table of Contents
- [Code of Conduct](#Code-of-Conduct)
- [How to Contribute](#How-to-Contribute)
- [Bug Reports](#Bug-Reports)
- [Feature Requests](#Feature-Requests)
- [Pull Requests](#Pull-Requests)
- [Running Tests](#Running-Tests)
- [Code Style](#Code-Style)
- [License](#License)
- [Community](#Community)
## Code of Conduct
By participating in this project, you agree to abide by our [Code of Conduct](https://github.com/lignum-vitae/goombay/blob/master/docs/CODE_OF_CONDUCT.md). Please take a moment to familiarize yourself with it.

## How to Contribute
### Bug Reports
If you find a bug or unexpected behavior, please open an issue in the GitHub repository. When reporting a bug, provide the following information:

- A description of the problem.
- Steps to reproduce the issue (if applicable).
- Any relevant error messages.
- The version of the library you're using.
- The Python version you're using.

### Feature Requests
If you have an idea for a new feature or enhancement, please open an issue describing the feature and why you think it would be useful. We encourage open discussions before starting to code a new feature.

### Pull Requests
To contribute code:

#### 1. Open a new Issue following the above-mentioned guidelines
#### 2. Fork the repository to your own GitHub account.
#### 3. Clone your fork locally:
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

#### 5. Checkout your branch for changes:
##### Choose ONE of the following commands
```nginx
# Creates a new branch that stays in sync with the develop branch of the main repository
git checkout -b feature-name upstream/develop

# Checks out existing branch if you already have a branch locally
git checkout feature-name
```
#### 6. Make your changes in your local repository and run in your local environment.
```nginx
# Downloads project as editable, which allows local imports
python -m pip install -e .

# Files may be run as a module
python -m goombay.algorithms.editdistance

# Or as a script
python editdistance.py
```
#### 7. Test your changes. Make sure all tests pass. (See Running Tests).
#### 8. Commit your changes with a descriptive commit message:
```nginx
# Gets latest changes from main biobase project if you've set up an upstream branch as detailed above
git fetch upstream
# We recommend individually adding each file with modifications
git add <filename>
# Commit files after all files with modifications have been added
git commit -m "Add feature: description of change"`
```
#### 🚨 Using git add .
While `git add .` is convenient for adding all modified files, it can lead to messy commits. Consider using it only when:
- You've reviewed all changes
- You're certain about each modification
- You've checked git status first
- Your .gitignore is properly configured
#### 9. Rebase Your Development Branch with the Latest Upstream changes
```nginx
# Make sure all is committed (or stashed) as necessary on this branch
git rebase -i upstream/develop feature-name
```
You may need to resolve conflicts that occur when both a file on the development trunk and one of the files in your branch have been changed.
Edit each conflicting file to resolve the differences, then continue the rebase.
Each file will need to be "added" to mark the conflict as resolved:
```nginx
# Resolve conflicts in each file, then:
git add <resolved-filename>
git rebase --continue
```
#### 10. Push your branch to your fork on GitHub:
```nginx
git push -f origin feature-name
```
#### 11. Open a Pull Request (PR) from your branch to the main branch of the original Goombay repository.
- You may need to click the `compare across forks` link under the `Compare changes` header that populates
  when you click `New pull request` to see your local repo fork.

#### 12. In your PR description, include:

- A summary of the changes.
- Any relevant issue numbers (e.g., fixes #123).
- Information about tests and validation.

## Running Tests
To ensure your changes work correctly, you can run the tests before submitting a PR.

Install dependencies (make sure you have numpy installed as well):

```nginx
pip install goombay
```
Run all tests using the following command from the root directory:

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
We follow standard Python conventions (PEP 8) for code style. Some additional notes:

- Use meaningful variable and function names.
- Keep lines of code under 80 characters where possible.
- Ensure all new features are covered by tests.
- Make sure to update documentation if your changes affect the usage or API.
## License
By contributing to Goombay, you agree that your contributions will be licensed under the MIT License, as outlined in the [LICENSE](https://github.com/lignum-vitae/goombay/blob/master/LICENSE) file.

## Community
We encourage contributions from everyone, and we strive to maintain a welcoming and inclusive community. If you have any questions, need help, or want to discuss ideas, feel free to reach out via issues or the repository discussions.

Thank you for contributing to Goombay! Your help improves the project for everyone!

