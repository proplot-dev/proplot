========
Tutorial
========

Below is a showcase of ProPlot's features. Once you have a general understanding of the API, see the :ref:`documentation` of ProPlot's individual classes and functions.

.. include:: showcase.rst

.. WARNING:
   I don't know how to write an nbconvert template so instead
   I used the following regex to delete cell text output:
   :%s/.. parsed-literal::\n\n.\+\_.\{-}\n\n//g
   And used the following to transform literal double backticks into
   sphinx links:
   :%s/``\(\~.\{-}\)``/`\1`/g
   :%s/:ref:``\(.\{-}\)``/:ref:`\1`/g
   Takes just a couple seconds.

