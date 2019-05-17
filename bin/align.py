"""
Align for sequence.

See more on https://en.wikipedia.org/wiki/Sequence_alignment.
"""


class CheckSequence(object):
    """Check sequence type."""

    sequencetype = set()

    def __init__(self, name):
        """Init name."""
        self.name = name

    def __get__(self, instance, cls):
        """Get value."""
        if instance is None:
            return self
        else:
            return instance.__dict__[self.name]

    def __set__(self, instance, value):
        """Set value."""
        value = value.strip().upper()
        if not set(value).issubset(self.sequencetype):
            raise TypeError("Sequence %s must be in sequencetype:%s!" %
                            (value, self.sequencetype))
        instance.__dict__[self.name] = value

    def __delete__(self, instance):
        """Delete value."""
        del instance.__dict__[self.name]


class ProteinSequence(CheckSequence):
    """Check whether is protein sequence."""

    sequencetype = set('ARNDCQEGHILKMFPSTWYVX')


class DNASequence(CheckSequence):
    """Check whether is DNA sequence."""

    sequencetype = set('ATCGN')


class Align(object):
    pass
