import os

def check_file_exists(file_path):
    """
    Check whether the file at file_path exists.

    :param file_path: full path to the file checked.
    :returns: a boolean indicating whether the file exists.
    """

    if os.path.isfile(file_path):
        return True
    else:
        return False
