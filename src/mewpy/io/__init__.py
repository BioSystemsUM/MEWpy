from pathlib import Path
from typing import Union, TYPE_CHECKING

from mewpy.simulation import get_container, get_simulator

from .director import Director
from .reader import Reader
from .writer import Writer

from .engines import (Engines,
                      BooleanRegulatoryCSV,
                      CoExpressionRegulatoryCSV,
                      TargetRegulatorRegulatoryCSV,
                      CobraModel,
                      ReframedModel,
                      JSON,
                      RegulatorySBML,
                      MetabolicSBML)


if TYPE_CHECKING:
    from io import TextIOWrapper

    from mewpy.germ.models import Model, RegulatoryModel, MetabolicModel
    from cobra import Model as Cobra_Model
    from reframed import CBModel as Reframed_Model


def load_sbml_container(filename, flavor='reframed'):
    if flavor == 'reframed':
        from reframed.io.sbml import load_cbmodel
        model = load_cbmodel(filename)
    elif flavor == 'cobra':
        from cobra.io import read_sbml_model
        model = read_sbml_model(filename)
    else:
        raise ValueError(f"{flavor} is not a recognized flavor")
    container = get_container(model)
    return container


def load_sbml_simulator(filename, flavor='reframed', envcond=None):
    if flavor == 'reframed':
        from reframed.io.sbml import load_cbmodel
        model = load_cbmodel(filename)
    elif flavor == 'cobra':
        from cobra.io import read_sbml_model
        model = read_sbml_model(filename)
    else:
        raise ValueError(f"{flavor} is not a recognized flavor")
    simul = get_simulator(model, envcond=envcond)
    return simul


def load_gecko_simulator(filename, flavor='reframed', envcond=None):
    if flavor == 'reframed':
        from mewpy.model.gecko import GeckoModel
        model = GeckoModel(filename)
    elif flavor == 'cobra':
        from geckopy.gecko import GeckoModel
        model = GeckoModel(filename)
    else:
        raise ValueError(f"{flavor} is not a recognized flavor")
    simul = get_simulator(model, envcond=envcond)
    return simul


# A common entry point for all models.
# The following "read" and "write" functions available in this module are just wrappers for a single reader/writer
# or multiple readers/writers using the director
def read_model(*readers: Reader,
               warnings: bool = True) -> Union['Model', 'RegulatoryModel', 'MetabolicModel']:
    """
    Reading a GERM model encoded into one or more file types (e.g. sbml, csv, cobrapy, reframed, json, etc).
    It can return a metabolic, regulatory or metabolic-regulatory model from multiple files.

    A reader must be provided for each file type.
    Reading will take place according to the order and settings of the each reader.

    All files are closed upon reading or failure

    :param readers: Multiple Reader instances that will be used to read multiple file types into a single GERM model
    :param warnings: Whether to launch warnings found during reading
    :return: mewpy metabolic, regulatory or both model
    """

    reader_director = Director(*readers)

    model = reader_director.read()

    if warnings:
        reader_director.warn()

    return model


def read_sbml(io: Union[str, Path, 'TextIOWrapper'],
              metabolic: bool = True,
              regulatory: bool = True,
              warnings: bool = True) -> Union['Model', 'RegulatoryModel', 'MetabolicModel']:
    """
    Reading a GERM model encoded into a SBML file.
    It can return a metabolic, regulatory or metabolic-regulatory model from the SBML file according to the metabolic
    and regulatory flags.

    Only one SBML file is accepted. Thus, if the SBML file does not encode an SBML-qual or SBML-fbc plugin,
    but either metabolic or regulatory reading are requested, an IO error will be raised, respectively.

    The SBML file is closed upon reading or failure

    :param io: A valid string path or IO.
    :param metabolic: Whether to read the metabolic share of the SBML file, namely the fbc plugin
    :param regulatory: Whether to read the regulatory share of the SBML file, namely the qual plugin
    :param warnings: Whether to launch warnings found during reading
    :return: mewpy metabolic, regulatory or both model
    """

    readers = []

    if metabolic:
        metabolic_reader = Reader(engine=MetabolicSBML,
                                  io=io)

        readers.append(metabolic_reader)

    if regulatory:
        regulatory_reader = Reader(engine=RegulatorySBML,
                                   io=io)

        readers.append(regulatory_reader)

    if not readers:
        raise OSError('Nothing to read')

    return read_model(*readers, warnings=warnings)


def read_csv(io: Union[str, Path, 'TextIOWrapper'],
             boolean: bool = True,
             co_expression: bool = False,
             target_regulator: bool = False,
             warnings: bool = True,
             **kwargs) -> Union['Model', 'RegulatoryModel']:
    """
    Reading a mewpy regulatory model encoded into a CSV file.
    It can only return a regulatory model from the CSV file.

    Only one CSV file is accepted. The CSV file type can be explicitly set using the boolean (default), co_expression
    or target_regulator flags. Consult the Engines enumerator in mewpy.io.engines for further detail on this CSV file
    types.

    The CSV file is closed upon reading or failure

    :param io: A valid string path or IO.
    :param boolean: Whether the file is a Boolean-based regulatory CSV file
    :param co_expression: Whether the file is a CoExpression (co-activating and co-repressing) regulatory CSV file
    :param target_regulator: Whether the file is a Target-Regulator interaction-based regulatory CSV file
    :param warnings: Whether to launch warnings found during reading
    :return: mewpy regulatory model
    """

    readers = []

    if boolean:
        boolean_reader = Reader(engine=BooleanRegulatoryCSV,
                                io=io,
                                **kwargs)

        readers.append(boolean_reader)

    if co_expression:
        regulatory_reader = Reader(engine=CoExpressionRegulatoryCSV,
                                   io=io,
                                   **kwargs)

        readers.append(regulatory_reader)

    if target_regulator:
        regulatory_reader = Reader(engine=TargetRegulatorRegulatoryCSV,
                                   io=io,
                                   **kwargs)

        readers.append(regulatory_reader)

    if not readers:
        raise OSError('Nothing to read')

    if len(readers) > 1:
        raise OSError('read_csv only accepts one file at a time')

    return read_model(*readers, warnings=warnings)


def read_cbmodel(io: Union['Cobra_Model', 'Reframed_Model'],
                 cobrapy: bool = True,
                 reframed: bool = False,
                 warnings: bool = True) -> Union['Model', 'MetabolicModel']:
    """
    Reading a mewpy metabolic model encoded into a Constraint-Based metabolic model from Cobrapy or Reframed.
    It can only return a metabolic model from the cobra model.

    Only one cobra model. The cobra model platform can be explicitly set using the cobrapy (default) or reframed flags.
    Consult the Engines enumerator in mewpy.io.engines for further detail on this cobra model types.

    :param io: A valid cobra model.
    :param cobrapy: Whether the cobra model is a cobrapy model
    :param reframed: Whether the cobra model is a reframed model
    :param warnings: Whether to launch warnings found during reading
    :return: mewpy metabolic model
    """

    readers = []

    if cobrapy:
        boolean_reader = Reader(engine=CobraModel,
                                io=io)

        readers.append(boolean_reader)

    if reframed:
        regulatory_reader = Reader(engine=ReframedModel,
                                   io=io)

        readers.append(regulatory_reader)

    if not readers:
        raise OSError('Nothing to read')

    if len(readers) > 1:
        raise OSError('read_cbmodel only accepts one cbmodel at a time')

    return read_model(*readers, warnings=warnings)


def read_json(io: Union[str, Path, 'TextIOWrapper'],
              warnings: bool = True) -> Union['Model', 'MetabolicModel', 'RegulatoryModel']:
    """
    Reading a GERM model encoded into a JSON file.
    It can return a metabolic, regulatory or metabolic-regulatory model from the JSON file according to the JSON file.

    Only one JSON file is accepted. The GERM model is built according to the JSON file content.

    The JSON file is closed upon reading or failure

    :param io: A valid string path or IO.
    :param warnings: Whether to launch warnings found during reading
    :return: mewpy metabolic, regulatory or both model
    """

    json_reader = Reader(engine=JSON,
                         io=io)

    return read_model(json_reader, warnings=warnings)


def write_model(*writers: Writer,
                warnings=True) -> Union['Model', 'RegulatoryModel', 'MetabolicModel']:
    """
    Writing a GERM model into one or more file types (e.g. sbml, csv, cobrapy, reframed, json, etc).
    It can write a metabolic, regulatory or metabolic-regulatory model to multiple files.

    A writer must be provided for each file type.
    Reading will take place according to the order and settings of the each writer.

    All files are closed upon writing or failure

    :param writers: Multiple Writer instances that will be used to write multiple file types from a single GERM model
    :param warnings: Whether to launch warnings found during reading
    :return: mewpy metabolic, regulatory or both model
    """

    writer_director = Director(*writers)

    model = writer_director.write()

    if warnings:
        writer_director.warn()

    return model
