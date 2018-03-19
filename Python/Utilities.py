from ConfigParser import SafeConfigParser
import logging.config


def config_dict(filename, tokenize=None):
    """
    A simple way to get a dictionary from a config file.
    If a list of strings is passed in tokenize, then those items will be split
    by commas and stripped of whitespace.
    """
    logging.config.fileConfig('logging.ini')

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
        logging.error("Could not find file %s", filename)
        raise

    if return_dict:
        for s, section in return_dict.iteritems():
            del section['__name__']
            if tokenize:
                assert isinstance(tokenize, list)
                for v, value in section.iteritems():
                    if ',' in value and v in tokenize:
                        section[v] = [i.strip() for i in value.split(',')]
        return return_dict
    else:
        logging.error("Config empty, something has gone wrong.")


def config_display(cfg_dict):
    return_str = "Listing config..."
    for s, section in cfg_dict.iteritems():
        return_str += "\n" + s + ":"
        for v, value in section.iteritems():
           return_str += "\n--" + v + " = " + str(value)

    return return_str

