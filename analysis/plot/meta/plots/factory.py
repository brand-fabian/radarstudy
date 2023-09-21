import typing
from abc import ABC, abstractmethod

from meta.statistics import StatisticsFactory


class AnalysisFactory:
    STATS_FACTORY: typing.Optional[StatisticsFactory] = None
    LANGUAGE: typing.Literal["en", "de"] = "en"

    @classmethod
    def instance(cls, name: str, *args, **kwargs) -> typing.Any:
        for subt in cls.__subclasses__():
            if subt.__name__ == name:
                val = subt.get(*args, **kwargs)
                return val
        return None

    @classmethod
    def types(cls) -> typing.List[str]:
        return list(map(lambda x: x.__name__, cls.__subclasses__()))

    @classmethod
    def create_stats_factory(cls):
        cls.STATS_FACTORY = StatisticsFactory()

    @classmethod
    @abstractmethod
    def get(cls, *args, **kwargs):
        pass
