from collections import defaultdict
from typing import Union, TYPE_CHECKING, Tuple

from mewpy.germ.models import Model
from .engines import JSON, MetabolicSBML

if TYPE_CHECKING:
    from mewpy.germ.models import Model, MetabolicModel, RegulatoryModel
    from .builder import Builder
    from .reader import Reader
    from .writer import Writer


class Director:

    def __init__(self, *builders: Union['Builder', 'Reader', 'Writer']):
        """
        Director is responsible for managing builders. It gives instructions to the builders on how to read or write
        a multi-type model properly

        To understand how reading and writing are proceeded, see the read and write methods

        :param builders: builders for reading models, files or IOs. builders for writing files using mewpy models
        """
        self._builders = builders

    @property
    def builders(self) -> Tuple[Union['Builder', 'Reader', 'Writer']]:
        """
        Returns the builders associated with this director
        :return: tuple of builders
        """
        return self._builders

    def read(self) -> Union['Model', 'RegulatoryModel', 'MetabolicModel']:
        """
        Reading a GERM model, namely metabolic, regulatory or both encoded into one or more file types.
        Reading is performed step-wise according to the builders order.
        :return: metabolic, regulatory or germ model
        """

        types = set()
        variables = defaultdict(set)

        for builder in self.builders:

            engine = builder.engine

            # First, opening the resource in read mode
            engine.open(mode='r')

            # Second, retrieving the model type
            types.add(engine.model_type)

            # Third, parsing the resource

            # parsing is responsible for two main operations:
            #   - reading a file or model with the adequate engine (pandas, python-lib-sbml, cobra, etc)
            #   and parse all important items to a data transfer object, the dto.

            #   - while reading the file or model, all variables and their types across the various files and models are
            #   collected and centralized under the engine.variables property
            engine.parse()

            # merging all variables and their types across the various files and models
            for variable, var_types in engine.variables.items():
                variables[variable].update(var_types)

            # Fourth, the resource opened by the engine are closed and can be used next
            engine.close()

        # Fifth, the model is created, as we now know the model types.
        model = Model.from_types(types=types, identifier='model')
        model_id = None
        model_name = None

        for builder in self.builders:

            engine = builder.engine

            if isinstance(engine, JSON):

                # json is a special case because of serialization

                model = engine.read()
                engine.clean()
                return model

            # Sixth, model building takes place.

            # The model is incremented with further containers and attributes, namely the stuff that each file/model
            # has to offer

            # In detail, when the read method of an engine receives a model object, it does not create a new one.
            # Instead, it adds new variables to the already available attributes and containers

            # In detail, when the read method of a engine receives a dictionary of a set of variables, it will create
            # the variables according to the variables types listed among all files.
            # That is, if var_01 has target, regulator and metabolite types registered across the various files/models,
            # the read method will create a TargetRegulatorMetabolite variable under a single instance
            # (a multi-type variable). This variable contains all attributes of a Target, Regulator and Metabolite
            engine.read(model=model, variables=variables)

            # update model id and name
            if isinstance(engine, MetabolicSBML):
                model_id = engine.dto.id
                model_name = engine.dto.name

            # Seventh, cleaning the data transfer object In detail, when the open and parse methods of an engine are
            # used, the data transfer object is populated with records. Then, it is no longer required for this
            # object to live in memory
            engine.clean()

        if model_id:
            model._id = model_id

        if model_name:
            model.name = model_name

        model.clean_history()
        return model

    def write(self):

        """
        Writing a GERM model, namely metabolic, regulatory or both to one or more file types.
        Writing is performed step-wise according to the builders order.
        :return:
        """
        for builder in self.builders:

            engine = builder.engine

            # First, opening the resource in write mode
            engine.open(mode='w')

            # Second, writing the model to the file
            engine.write()

            # Third, closing the resource
            engine.close()

            # Seventh, cleaning the data transfer object In detail, when the open and write methods of an engine are
            # used, the data transfer object is populated with records. Then, it is no longer required for this
            # object to live in memory
            engine.clean()

    def warn(self):
        """
        Launch warnings stored in the builders.
        :return:
        """
        # warnings and even small errors are all collected with partial pattern and then launched, as printing to the
        # console takes time

        for builder in self.builders:

            for warning in builder.warnings:
                warning()
