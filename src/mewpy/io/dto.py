from dataclasses import dataclass, field
from typing import Any, Dict, Union, TYPE_CHECKING, Tuple, Set

from pandas import DataFrame

try:
    # noinspection PyPackageRequirements
    from cobra import Model as Cobra_Model

except ImportError:
    Cobra_Model = str

try:
    # noinspection PyPackageRequirements
    from reframed import CBModel as Reframed_Model

except ImportError:

    Reframed_Model = str

from mewpy.mew.algebra import Symbolic, NoneAtom
from mewpy.mew.variables import Variable

if TYPE_CHECKING:
    from mewpy.mew.variables import Gene, Interaction, Metabolite, Reaction, Regulator, Target
    from mewpy.model import Model, MetabolicModel, RegulatoryModel


@dataclass
class History:

    """
    Data transfer object for history.
    History encoded into a sbml model
    """

    data: str = None
    creators: str = None


@dataclass
class FunctionTerm:

    """
    Data transfer object for function terms.
    Function term holds the symbolic expression of a given interaction between a target and a set of regulators.
    It also holds the resulting coefficient of this interaction
    """

    id: str = None
    symbolic: Symbolic = field(default_factory=NoneAtom)
    coefficient: int = 0


@dataclass
class CompartmentRecord:

    """
    Data transfer object for compartments.
    """

    id: str = 'e'
    name: str = 'external'


@dataclass
class VariableRecord:

    """
    Data transfer object for variables encoded into files.
    Variable record holds information regarding any multi-type variable, namely metabolic or regulatory.
    It can be extended to more types.
    It stores the many attributes that these variables can have encoded into multiple file types.
    """

    # ids, names, types, etc
    id: Any = field(default_factory=str)
    name: str = field(default_factory=str)
    aliases: set = field(default_factory=set)

    # other
    notes: Any = None
    annotation: Any = None

    # common attributes
    compartment: Any = None
    coefficients: set = field(default_factory=set)
    active_coefficient: float = 0.0
    constant: Any = None

    # metabolic attributes
    bounds: tuple = None
    charge: int = None
    formula: str = None
    genes: Dict[str, 'VariableRecord'] = field(default_factory=dict)
    gpr: FunctionTerm = field(default_factory=FunctionTerm)
    metabolites: Dict[str, 'VariableRecord'] = field(default_factory=dict)
    products: Dict[str, 'VariableRecord'] = field(default_factory=dict)
    reactants: Dict[str, 'VariableRecord'] = field(default_factory=dict)
    stoichiometry: Dict[str, Union[int, float]] = field(default_factory=dict)

    # regulatory attributes
    regulators: Dict[str, 'VariableRecord'] = field(default_factory=dict)
    target: 'VariableRecord' = None
    interactions: Dict[str, 'VariableRecord'] = field(default_factory=dict)
    function_terms: Dict[str, 'FunctionTerm'] = field(default_factory=dict)

    # TODO: it should be improved for further automation, so that one can just type to_variable(). For now,
    #  one still has to pass the correct attributes at each reading/model building
    def to_variable(self,
                    model: Union['Model', 'MetabolicModel', 'RegulatoryModel'],
                    types: Set[str],
                    **attributes) -> Tuple[Union['Gene', 'Interaction', 'Metabolite', 'Reaction', 'Regulator',
                                                 'Target'], str]:

        try:

            variable = model.get(self.id, None)

        except AttributeError:

            warning = f'{self.id} cannot be built properly. Generic Variable, build instead'

            return Variable(identifier=self.id), warning

        if variable is None:

            try:

                attributes['identifier'] = self.id

                variable = Variable.from_types(types=types, **attributes)

                return variable, ''

            except TypeError:

                warning = f'{self.id} cannot be built properly. Generic Variable, build instead'

                return Variable(identifier=self.id), warning

        else:

            try:

                variable.update(**attributes)

                return variable, ''

            except TypeError:

                warning = f'{self.id} cannot be built properly. Generic Variable, build instead'

                return variable, warning


@dataclass
class DataTransferObject:
    """
    Data transfer object for models encoded into files. This is the primary DTO associated with a builder.
    DTO record holds information regarding any multi-type model, namely metabolic or regulatory.
    It can be extended to more types.
    It stores the many attributes that these models can have encoded into multiple file types.
    """

    # cobra model
    cobra_model: Cobra_Model = None

    # reframed model
    reframed_model: Reframed_Model = None

    # CSV, TXT dataframe
    data_frame: DataFrame = field(default_factory=DataFrame)
    aliases_columns: list = field(default_factory=list)

    # SBML doc & model
    doc: Any = None
    model: Any = None

    # SBML plugins
    qual_plugin: Any = None
    fbc_plugin: Any = None

    # ids, names, versions, etc
    id: Any = None
    name: str = None
    level: Any = None
    version: Any = None
    history: Any = None

    # types
    types: set = None

    # common containers
    compartments: Dict[str, VariableRecord] = field(default_factory=dict)
    extracellular_reactions: Dict[str, VariableRecord] = field(default_factory=dict)
    variables: Dict[str, VariableRecord] = field(default_factory=dict)

    # metabolic containers
    extracellular_metabolites: Dict[str, VariableRecord] = field(default_factory=dict)
    genes: Dict[str, VariableRecord] = field(default_factory=dict)
    metabolites: Dict[str, VariableRecord] = field(default_factory=dict)
    objective: Dict[str, VariableRecord] = field(default_factory=dict)
    reactions: Dict[str, VariableRecord] = field(default_factory=dict)

    # regulatory containers
    interactions: Dict[str, VariableRecord] = field(default_factory=dict)
    regulators: Dict[str, VariableRecord] = field(default_factory=dict)
    targets: Dict[str, VariableRecord] = field(default_factory=dict)
