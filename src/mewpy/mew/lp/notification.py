from typing import Union, List, Tuple, Set, Dict, TYPE_CHECKING

if TYPE_CHECKING:
    from mewpy.mew.variables.coefficient import Coefficient


class Notification:

    def __init__(self,
                 content: Union[Dict, List, Tuple, Set, 'Coefficient'],
                 content_type: str,
                 action: str):
        """
        A notification object. It is used to notify the user about changes in the model.
        A notification contains a message and a payload. The message carries the content of the changes and the payload
        carries the information about the changes.

        :param content: the content of the notification aka the message for the linear problem to read
        :param content_type: the type of the content. It can be a variable, a reaction, a metabolite, etc. aka payload
        :param action: the action that triggered the notification. It can be an addition, a removal, etc. aka payload
        """
        content_types = ('reactions', 'metabolites', 'gprs', 'interactions', 'coefficients', 'objectives')
        actions = ('add', 'remove', 'set')

        if content_type not in content_types:
            raise ValueError(f'Wrong content_type, select one of the following: {content_types}')

        if action not in actions:
            raise ValueError(f'Wrong action, select one of the following: {actions}')

        self.content = content
        self.content_type = content_type
        self.action = action
