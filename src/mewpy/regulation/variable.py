
class RegulatoryVariable:

    def __init__(self, id=None, name=None, aliases=None):

        self._id = id
        self._name = name
        self._aliases = aliases

    def __str__(self):
        return self._id

    def __repr__(self):

        return '{}({}, {}, {}) at {}'.format(self.__class__.__name__, self._id, self._name, self._aliases, id(self))

    @property
    def id(self):
        return getattr(self, '_id', None)

    @property
    def name(self):
        return getattr(self, '_name', None)

    @property
    def aliases(self):
        return getattr(self, '_aliases', None)

    @id.setter
    def id(self, value):

        if value == self._id:
            pass

        elif not isinstance(value, str):
            raise TypeError("The ID must be a string. Use str() to convert")

        else:
            self._id = value

    @name.setter
    def name(self, value):

        if value == self._name:
            pass

        elif not isinstance(value, str):
            raise TypeError("The name must be a string. Use str() to convert")

        else:
            self._name = value

    @aliases.setter
    def aliases(self, value):

        if value == self._aliases:
            pass

        elif not isinstance(value, list):
            raise TypeError("The aliases must be a list. Use list() to convert")

        else:
            self._aliases = value
