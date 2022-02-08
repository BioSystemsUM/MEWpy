import random
import string

from reframed.core.cbmodel import CBModel
from reframed.external.cobrapy import to_cobrapy


def escher_maps():
    try:
        import escher
    except ImportError:
        raise RuntimeError("Escher is not installed.")

    maps = escher.list_available_maps()
    return [entry['map_name'] for entry in maps]


def randomString(stringLength=10):
    """Generate a random string of fixed length """
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))


def to_json(model, filename=None):
    import cobra
    if not filename:
        filename = randomString() + ".json"
    if isinstance(model, cobra.core.model.Model):
        c_model = model
    elif isinstance(model, CBModel):
        c_model = to_cobrapy(model)
    else:
        raise Exception
    cobra.io.save_json_model(c_model, filename)
    return filename


def remove_prefix(text, prefix='R_'):
    if text.startswith(prefix):
        return text[len(prefix):]
    return text


def build_escher(model=None, fluxes=None, fmt_func=remove_prefix, **kwargs):
    try:
        import escher
    except ImportError:
        raise RuntimeError("Escher is not installed.")

    js = None
    if model is None:
        map_name = 'e_coli_core.Core metabolism'
    elif isinstance(model, str) and model in escher_maps():
        map_name = model
    else:
        try:
            js = to_json(model)
        except Exception:
            map_name = 'e_coli_core.Core metabolism'

    if fluxes and fmt_func:
        fluxes = {fmt_func(r_id): val for r_id, val in fluxes.items()}
    else:
        fluxes = None
    if js:
        return escher.Builder(model_json=js, reaction_data=fluxes, **kwargs)
    else:
        return escher.Builder(map_name=map_name, reaction_data=fluxes, **kwargs)
