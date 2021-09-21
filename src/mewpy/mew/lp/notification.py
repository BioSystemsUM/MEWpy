from typing import Union, List, Tuple, Set, Dict, TYPE_CHECKING

if TYPE_CHECKING:
    from mewpy.mew.variables.coefficient import Coefficient


class Notification:

    def __init__(self,
                 content: Union[Dict, List, Tuple, Set, 'Coefficient'],
                 content_type: str,
                 action: str):

        content_types = ('reactions', 'metabolites', 'gprs', 'interactions', 'coefficients', 'objectives')
        actions = ('add', 'remove', 'set')

        if content_type not in content_types:
            raise ValueError(f'Wrong content_type, select one of the following: {content_types}')

        if action not in actions:
            raise ValueError(f'Wrong action, select one of the following: {actions}')

        self.content = content
        self.content_type = content_type
        self.action = action
