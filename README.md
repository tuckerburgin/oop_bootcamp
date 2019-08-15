# oop_bootcamp

oop_examples.py contains several examples of basic object-oriented programming practices

---

The other files represent three iteratively more advanced approaches to using multiple different MD analysis packages
for the same task, handling their different syntaxes:

---

First, scripy.py handles the problem entirely manually.

Second, interface.py and interfaced_script.py adds an interfacing class to greatly simply the problem.

Third, interface.py, factory.py, and factory_script.py makes the problem arbitrarily extensible by allowing users
to define their own interface scripts without needing to modify the 'user' script (factory_script.py) (instead, they
would modify interface.py).

Finally, registry_interface.py, factory_registry.py, and registry_script.py abstracts the definition of an adapter from
the interface module, so that to add a new adapter all a user would have to do is make a new file for that particular
toolkit and add it to the interfaces/ directory (this functionality is broken but close).

---

All are called as: python [script] [tool], where [tool] is mdanalysis or mdtraj (or others, if you add them!)
