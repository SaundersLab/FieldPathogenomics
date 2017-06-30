Protocols
===========

Generating a new callset
-------------------------

**1. Generate a text file containing the library names**

    .. code-block:: bash

        echo 'LIB22849
              LIB23170
                 .
                 .
                 .
              LIB23262' > my_callset.txt

**2. Run the Callset script*

Either run this in the production area or you local dev area

     .. code-block:: bash

        ./scripts/Callset.sh my_callset.txt 10

The Callset script takes two parameters, the list of libraries to include in the callset (and it uses the file name as the callset name) and an integer number of parallel blocks to spilt the callset into.


Generating a new tree
-----------------------

Exactly the same as with Callset,

     .. code-block:: bash

        ./scripts/Tree.sh my_callset.txt
