from abc import ABC, abstractmethod
import numpy as np

class Calculator(ABC):
    def __init__(self) -> None:
        pass

    @abstractmethod
    def calculate(self):
        pass