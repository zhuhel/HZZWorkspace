from ConfigParser import SafeConfigParser
import logging
import os


def config_dict(filename, tokenize=None):
    """
    A simple way to get a dictionary from a config file.
    If a list of strings is passed in tokenize, then those items will be split
    by commas and stripped of whitespace.
    """
    config = SafeConfigParser()
    config.optionxform = str
    return_dict = None

    logging.info("Attempting to read file %s...", filename)
    try:
        with open(filename) as f:
            config.readfp(f)
            return_dict = config._sections
            logging.info("Success!")
    except IOError:
        logging.error("Could not find file %s. Consider using absolute paths in your config file.", filename)
        raise

    if return_dict:
        for s, section in return_dict.iteritems():
            del section['__name__']
            if tokenize:
                assert isinstance(tokenize, list)
                for v, value in section.iteritems():
                    if v in tokenize:
                        section[v] = [i.strip() for i in value.split(',')]
        return return_dict
    else:
        logging.error("Config empty, something has gone wrong.")


def config_display(cfg_dict):
    """
    A nice way to display the given config dictionary.
    :param cfg_dict: the input config dictionary, as in the config_dict function above
    :return: a formatted string with the sections and values
    """
    return_str = "Listing config..."
    for s, section in cfg_dict.iteritems():
        return_str += "\n" + s + ":"
        for v, value in section.iteritems():
            return_str += "\n-- " + v + " = " + str(value)
    return return_str


def check_np(np_dict):
    """
    Checks to see if the NP config file is appropriately designed.
    :param np_dict: The config dict coming from the config_dict function above
    :return: a boolean, True if the dict is fine, False if not
    """
    for s, section in np_dict.iteritems():
        var_sum = 0
        for v, value in section.iteritems():
            if  "JET_JER_SINGLE_NP" in v:
                continue
            if "1up" in v:
                var_sum += 1
            elif "1down" in v:
                var_sum -= 1

        if var_sum == 0:
            return True
        else:
            return False