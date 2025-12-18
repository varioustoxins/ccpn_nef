"""
STAR File Representation Classes

This module defines the core classes for representing STAR/NEF file structures,
including both the current constructors and proposed enhanced constructors
that support declarative initialization.
"""

from collections import OrderedDict


class UnquotedValue(str):
    """String subclass for values that should not be quoted in STAR format."""
    pass


class NamedOrderedDict(OrderedDict):
    """Base class for named ordered dictionary structures."""

    def __init__(self, name=None):
        """
        Initialize a named ordered dictionary.
        
        Args:
            name (str, optional): Name identifier for this container
        """
        super(NamedOrderedDict, self).__init__()
        self.name = name

    def __str__(self):
        return '%s(name=%s)' % (self.__class__.__name__, self.name)

    def __repr__(self):
        return '%s(%s, name=%s)' % (
            self.__class__.__name__, 
            list(tt for tt in self.items()), 
            self.name
        )

    def addItem(self, tag, value):
        """
        Add an item to the container.
        
        Args:
            tag (str): Key identifier
            value: Value to store
            
        Raises:
            ValueError: If tag already exists
        """
        if tag in self:
            raise ValueError("%s: duplicate key name %s" % (self, tag))
        else:
            self[tag] = value


class StarContainer(NamedOrderedDict):
    """Base class for STAR containers (DataBlocks and SaveFrames)."""
    
    def multiColumnValues(self, columns):
        """
        Get tuple of orderedDict of values for columns.
        
        Args:
            columns: Column names to retrieve
            
        Returns:
            Tuple or OrderedDict of values, or None if no matches
        """
        valueDict = OrderedDict((x, self.get(x)) for x in columns)
        testSet = set(x for x in valueDict.values() if x is not None)
        if not testSet:
            return None
        # Implementation details omitted for brevity
        return valueDict


class Loop:
    """
    Loop for general STAR object tree.
    
    Attributes:
        name (str): Loop identifier
        columns (list): List of string column headers  
        data (list): List of rows, where rows are OrderedDicts
    """

    def __init__(self, name=None, columns=None, data=None):
        """
        Initialize a Loop.
        
        Current constructor:
            Loop(name=None, columns=None)
            
        Enhanced constructor:
            Loop(name=None, columns=None, data=None)
        
        Args:
            name (str, optional): Loop identifier
            columns (list, optional): Column header names
            data (list, optional): List of row dictionaries (enhanced)
        """
        self.name = name
        
        if columns:
            self._columns = list(columns)
        else:
            self._columns = []
            
        # Enhanced: Support data initialization
        if data:
            self.data = list(data)
        else:
            self.data = []

    @property
    def columns(self):
        """Get column names as tuple."""
        return tuple(self._columns)

    def addRow(self, **kwargs):
        """
        Add a row to the loop.
        
        Args:
            **kwargs: Column values as keyword arguments
        """
        row = OrderedDict()
        for col in self.columns:
            row[col] = kwargs.get(col, None)
        self.data.append(row)

    def __str__(self):
        return '<%s:%s>' % (self.__class__.__name__, self.name)


class DataBlock(StarContainer):
    """
    DataBlock for general STAR object tree.
    
    Top-level container for STAR file data.
    """
    
    # Current constructor signature: DataBlock inherits NamedOrderedDict.__init__(name=None)
    # Enhanced constructor would be: DataBlock(name=None, children=None)

    def __init__(self, name=None, children=None):
        """
        Initialize a DataBlock.
        
        Current usage:
            block = DataBlock(name='data_name')
            block.addItem('key', value)
            
        Enhanced usage:
            block = DataBlock(name='data_name', children={
                'save_frame_1': SaveFrame(...),
                'key': 'value'
            })
        
        Args:
            name (str, optional): DataBlock name
            children (dict, optional): Initial child items (enhanced)
        """
        super(DataBlock, self).__init__(name=name)
        
        # Enhanced: Support children initialization
        if children:
            for key, value in children.items():
                self.addItem(key, value)

    tagPrefix = None  # Can be set in subclass instances


class SaveFrame(StarContainer):
    """
    SaveFrame for general STAR object tree.
    
    Container for related data items and loops.
    """
    
    def __init__(self, name=None, children=None):
        """
        Initialize a SaveFrame.
        
        Current usage:
            frame = SaveFrame(name='frame_name') 
            frame.addItem('key', value)
            
        Enhanced usage:
            frame = SaveFrame(name='frame_name', children={
                'sf_category': 'category_name',
                'loop_name': Loop(...)
            })
        
        Args:
            name (str, optional): SaveFrame name
            children (dict, optional): Initial child items (enhanced)
        """
        super(SaveFrame, self).__init__(name=name)
        
        # Enhanced: Support children initialization
        if children:
            for key, value in children.items():
                self.addItem(key, value)


class NmrDataBlock(DataBlock):
    """DataBlock for NMRSTAR/NEF object tree."""

    def __init__(self, name=None, children=None):
        """
        Initialize an NMR DataBlock.
        
        Enhanced constructor supporting children parameter.
        
        Args:
            name (str, optional): DataBlock name
            children (dict, optional): Initial child items
        """
        super(NmrDataBlock, self).__init__(name=name, children=children)

    def newSaveFrame(self, name, category):
        """
        Create and add a new NmrSaveFrame.
        
        Args:
            name (str): SaveFrame name
            category (str): SaveFrame category
            
        Returns:
            NmrSaveFrame: The created SaveFrame
        """
        # Implementation would add sf_category and sf_framecode automatically
        saveFrame = NmrSaveFrame(name, category=category)
        self.addItem(name, saveFrame)
        saveFrame.addItem('sf_category', category)
        saveFrame.addItem('sf_framecode', name)
        return saveFrame


class NmrSaveFrame(SaveFrame):
    """SaveFrame for NMRSTAR/NEF object tree."""

    def __init__(self, name=None, category=None, children=None):
        """
        Initialize an NMR SaveFrame.
        
        Current usage:
            frame = NmrSaveFrame(name='frame_name', category='category')
            
        Enhanced usage:
            frame = NmrSaveFrame(
                name='frame_name', 
                category='category',
                children={
                    'sf_category': 'category',
                    'data_item': 'value'
                }
            )
        
        Args:
            name (str, optional): SaveFrame name
            category (str, optional): SaveFrame category
            children (dict, optional): Initial child items (enhanced)
        """
        super(NmrSaveFrame, self).__init__(name=name, children=children)
        self.category = category

    @property
    def tagPrefix(self):
        """Prefix to use before item tags on output."""
        return '_%s.' % self.category if self.category else None

    def newLoop(self, name, columns):
        """
        Create and add a new Loop.
        
        Args:
            name (str): Loop name
            columns (list): Column names
            
        Returns:
            NmrLoop: The created Loop
        """
        loop = NmrLoop(name, columns)
        self.addItem(name, loop)
        return loop


class NmrLoop(Loop):
    """Loop for NMRSTAR/NEF object tree."""

    def __init__(self, name=None, columns=None, data=None):
        """
        Initialize an NMR Loop.
        
        Enhanced constructor supporting data parameter.
        
        Args:
            name (str, optional): Loop name
            columns (list, optional): Column names
            data (list, optional): Initial row data
        """
        super(NmrLoop, self).__init__(name=name, columns=columns, data=data)


# Type aliases for convenience
DataBlock = NmrDataBlock
SaveFrame = NmrSaveFrame
Loop = NmrLoop


# Example usage demonstrating both approaches:
def example_current_usage():
    """Example using current constructor pattern."""
    
    # Create empty structures
    block = NmrDataBlock(name='nef_my_project')
    
    # Add items incrementally
    frame = NmrSaveFrame(name='nef_meta_data', category='nef_nmr_meta_data')
    frame.addItem('sf_category', 'nef_nmr_meta_data')
    frame.addItem('format_name', 'nmr_exchange_format')
    
    # Add loop
    loop = NmrLoop(name='nef_sequence', columns=['index', 'chain_code', 'residue_name'])
    loop.addRow(index=1, chain_code='A', residue_name='ALA')
    loop.addRow(index=2, chain_code='A', residue_name='VAL')
    
    frame.addItem('nef_sequence', loop)
    block.addItem('nef_meta_data', frame)
    
    return block


def example_enhanced_usage():
    """Example using enhanced constructor pattern."""
    
    # Create complete structure declaratively
    block = NmrDataBlock(
        name='nef_my_project',
        children={
            'nef_meta_data': NmrSaveFrame(
                name='nef_meta_data',
                category='nef_nmr_meta_data',
                children={
                    'sf_category': 'nef_nmr_meta_data',
                    'format_name': 'nmr_exchange_format',
                    'nef_sequence': NmrLoop(
                        name='nef_sequence',
                        columns=['index', 'chain_code', 'residue_name'],
                        data=[
                            {'index': 1, 'chain_code': 'A', 'residue_name': 'ALA'},
                            {'index': 2, 'chain_code': 'A', 'residue_name': 'VAL'}
                        ]
                    )
                }
            )
        }
    )
    
    return block


if __name__ == "__main__":
    # Test both approaches
    current_block = example_current_usage()
    enhanced_block = example_enhanced_usage()
    
    print("Current approach:")
    print(current_block)
    print("\nEnhanced approach:")
    print(enhanced_block)