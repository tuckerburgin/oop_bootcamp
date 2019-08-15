"""
Examples of OOP practices in Python 3
"""

import abc      # abstract base class

# ----------------------------------------------------------------------------------- #
#                                       Classes                                       #
# ----------------------------------------------------------------------------------- #

# class Item(object):   # For python 2 or 3
class Item:             # For python 3 only
    _count = 10         # class attribute; shared among ALL instances of this class (usually avoid)
    def __init__(self, name=''):                    # __init__ method
        self.name = name                            # instance attribute (or just "attribute")
        self._private_name = 'don\'t access'        # "private" attribute (solely conventional, says "this is internal")
        self.__really_private = 'can\'t access'     # actually cannot be accessed (directly)
        Item._count += 1                            # increment _count for ALL Item objects
        #self._count += 1                           # increment _count for only the present Item object's count

    @classmethod        # means we should change 'self' to 'cls' (convention); this is shared among all instances (rare)
    def print_count(cls):
        print('Count =', cls._count)

    @staticmethod       # means stricly a subset of this class (rare)
    def print_date():   # does not take 'self'
        print('Today is ...')
        print('I don\'t know my class or instance')

    def calculate_price(self):
        print('Calculate price in Item')


class MyAbstract(abc.ABC):
    @abc.abstractmethod     # means that this method must be defined in all subclasses that inherits from this class
    def my_abstract_method(self):
        print('Abstract method in a superclass')    # code is unreachable except by super() call in a subclass


class Cart:                 # example of aggregation (Cart 'has' Item objects)
    def __init__(self):
        self.items = []     # list of Item objects 'in' Cart

# ----------------------------------------------------------------------------------- #
#                                      Subclasses                                     #
# ----------------------------------------------------------------------------------- #

class BookItem(Item):
    def __init__(self, name, author):
        self.author = author
        super().__init__(name)                  # calls parent class (Item) (python 3 only)
        #super(BookItem, self).__init__(name)   # calls parent class (Item) (python 2 or 3)

    def get_author(self):
        print(self.author)

    def calculate_price(self):                  # overrides method with same name in parent class Item
        print('Calculate price in BookItem')


class SubAbstract(MyAbstract):
    def my_abstract_method(self):
        print('Abstract method in a subclass')
        super().my_abstract_method()    # calls my_abstract_method from superclass

# ----------------------------------------------------------------------------------- #

book_item = Item('Cats')

print(book_item.name)
print(book_item._private_name)
# print(book_item.__really_private)     # throws AttributeError

book_item.print_count()     # prints 11

another_item = Item('Dogs')
another_item.print_count()  # prints 12; Item._count was already incremented one by book_item.__init__()

Item.print_count()          # class method without specific copy of the class; prints 12 because __init__ not called

my_book = BookItem('Frogs', 'Sam Ellis')
my_book.print_count()       # prints 13 even though BookItem doesn't directly have this method (inherits it)
my_book.get_author()
print(my_book.name)         # only works because name attribute was defined in super() call in BookItem def __init__
my_book.calculate_price()   # uses BookItem version of this function

#an_abc = MyAbstract()      # throws TypeError; cannot instantiate abstract classes with an abstractmethod
an_abc = SubAbstract()
an_abc.my_abstract_method()

for i in (book_item, my_book):      # example of polymorphism
    i.calculate_price()             # book_item and my_book are of different type, but both have calculate_price
