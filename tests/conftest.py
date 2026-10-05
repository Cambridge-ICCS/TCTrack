"""Pytest configuration."""


def pytest_ignore_collect(collection_path):
    """Ignore the integration tests unless they are specifically called."""
    return collection_path.parts[-1] == "integration"
