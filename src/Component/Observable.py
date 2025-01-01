class Observable:
    def __init__(self):
        self._observers = []

    def addObserver(self, observer):
        """Add an observer."""
        if observer not in self._observers:self._observers.append(observer)

    def removeObserver(self, observer):
        """Remove an observer."""
        self._observers.remove(observer)
    def notifyObservers(self,**kwargs):
        for observer in self._observers: observer.update(**kwargs)



