from typing import Any, TYPE_CHECKING


if TYPE_CHECKING:
    from mewpy.mew.variables import Variable


def coefficients_setter(instance: 'Variable', value: Any):
    if hasattr(instance, '_bounds'):

        # if it is a reaction, bounds must be set
        instance._bounds = tuple(value)

    # if it is a metabolite, the bounds coefficient of the exchange reaction must be returned
    elif hasattr(instance, 'exchange_reaction'):

        if hasattr(instance.exchange_reaction, '_bounds'):
            instance.exchange_reaction._bounds = tuple(value)

    else:
        instance._coefficients = tuple(value)

    if instance.model:
        instance.model.notify()
