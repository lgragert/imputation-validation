# Release Guide (PyPI)

## Versioning

Before releasing, update all of these to the same version:

- `pyproject.toml` (`[project].version`)
- `setup.py` (`version=`)
- `voihla/__init__.py` (`__version__`)

Current planned release: `0.1.2`.

## One-Time Setup

Install build and upload tools in your active environment:

```bash
python -m pip install --upgrade build twine
```

Create a PyPI API token:

1. In PyPI account settings, create a token with scope for this project.
2. Save the token somewhere safe.

Optional: create `~/.pypirc`:

```ini
[pypi]
  username = __token__
  password = pypi-<your-token>

[testpypi]
  repository = https://test.pypi.org/legacy/
  username = __token__
  password = pypi-<your-testpypi-token>
```

## Release Steps

From repository root:

1. Clean old artifacts:

```bash
rm -rf build dist && find . -maxdepth 1 -name '*.egg-info' -exec rm -rf {} +
```

2. Build sdist + wheel:

```bash
python -m build
```

3. Validate package metadata:

```bash
python -m twine check dist/*
```

4. Upload to TestPyPI first (recommended):

```bash
python -m twine upload --repository testpypi dist/*
```

5. Verify install from TestPyPI:

```bash
python -m pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple voihla==0.1.2
```

6. Upload to PyPI:

```bash
python -m twine upload dist/*
```

## Notes

- PyPI does not allow overwriting an existing version.
- If upload fails due to version conflict, bump version and rebuild.
- Prefer tagging releases in git, for example: `v0.1.2`.
