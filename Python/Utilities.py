from ConfigParser import SafeConfigParser

def configDict(filename, tokenize=None):
    """
    A simple way to get a dictionary from a config file.
    If a list of strings is passed in tokenize, then those items will be split
    by commas and stripped of whitespace.
    """
    config = SafeConfigParser()
    config.optionxform = str
    with open(filename) as f:
        config.readfp(f)
        config_dict = config._sections
        
    for s, section in config_dict.iteritems():
        del section['__name__']
        if tokenize:
            assert isinstance(tokenize, list)
            for v, value in section.iteritems():
                if ',' in value and v in tokenize:
                    section[v] = [i.strip() for i in value.split(',')]
    return config_dict
