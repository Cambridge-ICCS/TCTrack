"""Temporary placeholder for TCTrack code as we set up repository boilerplate."""


def get_message(name: str = "my repo") -> str:
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
    hello_message = get_message("TCTrack")
    print(hello_message)
