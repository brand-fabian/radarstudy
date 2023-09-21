import typing

from .Serializers import FILE_TYPE, Serializer


class SerializerFactoryException(Exception):
    pass


def _get_types(
    root: typing.Type[Serializer],
) -> typing.Dict[FILE_TYPE, typing.Type[Serializer]]:
    ret_val = {}
    if root.provides() is not None:
        ret_val = {root.provides(): root}

    for subt in root.__subclasses__():
        ret_val = {**ret_val, **_get_types(subt)}
    return ret_val


class SerializerFactory:
    @classmethod
    def make_provider(cls, functionality_type: FILE_TYPE, *args, **kwargs):
        providers = _get_types(Serializer)
        if functionality_type in providers:
            return providers[functionality_type](*args, **kwargs)

        raise SerializerFactoryException(
            "missing provider for functionality {}".format(functionality_type)
        )

    @classmethod
    def get_available(cls):
        return list(_get_types(Serializer).keys())
