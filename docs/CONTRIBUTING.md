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

A description of the problem.
Steps to reproduce the issue (if applicable).
Any relevant error messages.
The version of the library you're using.
The Python version you're using.

### Feature Requests
If you have an idea for a new feature or enhancement, please open an issue describing the feature and why you think it would be useful. We encourage open discussions before starting to code a new feature.

### Pull Requests
To contribute code:

1. **Fork the repository** to your own GitHub account.
2. **Clone your fork** locally:
```nginx
git clone https://github.com/YOUR_USERNAME/goombay.git
```
3. **Create a new branch** for your changes:
```nginx
git checkout -b feature-name
```
4. **Make your changes** in your local repository.
5. **Test your changes**. Make sure all tests pass. (See Running Tests).
6. **Commit your changes** with a descriptive commit message:
```nginx
git commit -m "Add feature: description of change"
```
7. **Push your branch** to your fork on GitHub:
```nginx
git push origin feature-name
```
8. **Open a Pull Request (PR)** from your branch to the main branch of the original Goombay repository.

9. **In your PR description**, include:

- A summary of the changes.
- Any relevant issue numbers (e.g., fixes #123).
- Information about tests and validation.

## Running Tests
To ensure your changes work correctly, you can run the tests before submitting a PR.

Install dependencies (make sure you have numpy installed as well):

```nginx
pip install goombay
```
Run all tests using the following command:

```nginx
python -m unittest discover tests
```
or if you're using py:
```nginx
py -m unittest discover tests
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

