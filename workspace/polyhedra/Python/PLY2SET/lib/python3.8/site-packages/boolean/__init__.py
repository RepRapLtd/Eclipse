import re

__version__ = '1.1.0'


def boolean(value):
    """ Main function. Receives string and return boolean values """

    if re.match('^True|true|t|yes|y|1$', value):
        return True

    if re.match('^False|false|f|no|n|0$', value):
        return False

    return None
