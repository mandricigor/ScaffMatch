'''
Created on Jul 25, 2014

@author: igor
'''



class Settings(object):
    _instance = None
    
    _settings = {
                "aligner": "bowtie2",
                "bundle_size": 1,
                "key_size": 2
                }
    
    
    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(Settings, cls).__new__(
                                cls, *args, **kwargs)
        return cls._instance
    
    def update(self, new_settings):
        self._settings.update(new_settings)
        
    def get(self, key):
        return self._settings.get(key)
    
    
    def set(self, key, value):
        self._settings[key] = value
           
    def __str__(self):
        return "\n".join(["%s: %s" % x for x in self._settings.iteritems()])
