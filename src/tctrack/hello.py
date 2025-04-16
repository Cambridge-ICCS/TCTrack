"""Temporary placeholder for TCTrack code as we set up repository boilerplate."""


def hello(name: str = "my repo") -> str:
    """
    Return a hello statement with the repository name if supplied.

    Parameters
    ----------
    name: str | None
        The name of the repository supplied as a string.

    Returns
    -------
    str
        A hello message.
    """
    return f"Hello from {name}!"


if __name__ == "__main__":
    message = hello("TCTrack")
    print(message)
