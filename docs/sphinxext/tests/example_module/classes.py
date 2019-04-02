__all__ = ['Spam', 'Egg']


class BaseSpam(object):
    """
    Base class for Spam
    """

    def eat(self, time):
        """
        Eat some spam in the required time.
        """
        pass

    def buy(self, price):
        """
        Buy some MOAR spam.
        """
        pass


class Spam(BaseSpam):
    """
    The main spam
    """
    pass


class Egg(object):
    """
    An egg (no inheritance)
    """

    def eat(self, time):
        """
        Eat some egg in the required time.
        """
        pass

    def buy(self, price):
        """
        Buy some MOAR egg.
        """
        pass

    @property
    def weight(self):
        """
        The weight of an egg
        """
        return 0
