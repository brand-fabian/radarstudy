import functools
import logging
import os
import typing

from .serializers import FILE_TYPE, SerializerFactory

logger = logging.getLogger(__name__)


class CheckpointException(Exception):
    pass


class NoSerializerAvailable(CheckpointException):
    pass


class PathNotFound(CheckpointException):
    pass


class Checkpoint:

    """Enable or disable using the checkpoint functionality."""

    _is_enabled = True

    @classmethod
    def checkpoint(
        cls,
        *,
        prefix: typing.Optional[str] = None,
        name: typing.Optional[str] = None,
        output_type=FILE_TYPE.PICKLE
    ):
        """Checkpoint a function invocation by checking the filesystem.

        This decorator wraps an invocation of the passed function and will only
        call the underlying function if the data has not been generated yet. Otherwise,
        the returned data will be loaded from the filesystem.

        Parameters
        ----------
        prefix : str
                Output file path prefix (i.e. directory).
        name : str
            Output file name.
        output_type : FILE_TYPE
                    Type of the function return value and serializer used in the
                    application of the wrapper.

        Returns
        -------
        A Wrapper function that will write generated output from the function to a file
        with the given serializer, or _only_ load the file contents if the files already
        exist.
        """
        if prefix is None or name is None:
            raise PathNotFound("no path prefix or filename supplied.")

        if output_type not in SerializerFactory.get_available():
            raise NoSerializerAvailable(
                "no serializer exists for {}".format(output_type)
            )

        def decorator(fn):
            @functools.wraps(fn)
            def wrapper(*args, **kwargs):
                file_path = os.path.abspath(os.path.join(prefix, name))
                if name is None or not Checkpoint._is_enabled:
                    return fn(*args, **kwargs)
                else:
                    serializer = SerializerFactory.make_provider(output_type)
                    if os.path.exists(file_path):
                        logger.info("Resuming from checkpoint: {}".format(file_path))
                        return serializer.read(file_path)
                    else:
                        logger.info(
                            "Generating data for checkpoint: {}".format(file_path)
                        )
                        fn_out = fn(*args, **kwargs)
                        serializer.write(fn_out, file_path)
                        return fn_out

            return wrapper

        return decorator

    @classmethod
    def disable(cls):
        Checkpoint._is_enabled = False

    @classmethod
    def enable(cls):
        Checkpoint._is_enabled = True

    @classmethod
    def is_enabled(cls):
        return Checkpoint._is_enabled


class NamedCheckpoint(Checkpoint):

    _checkpoints = {}

    @classmethod
    def checkpoint(
        cls,
        *,
        prefix: typing.Optional[str] = None,
        name: typing.Optional[str] = None,
        output_type=FILE_TYPE.PICKLE
    ):

        NamedCheckpoint._checkpoints[name] = super().checkpoint(
            prefix=prefix, name=name, output_type=output_type
        )
        return NamedCheckpoint._checkpoints[name]

    @classmethod
    def keys(cls):
        return NamedCheckpoint._checkpoints.keys()

    @classmethod
    def items(cls):
        return NamedCheckpoint._checkpoints.items()

    @classmethod
    def values(cls):
        return NamedCheckpoint._checkpoints.values()
