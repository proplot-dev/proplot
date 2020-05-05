{% if referencefile %}
.. include:: {{ referencefile }}
{% endif %}

{{ objname }}
{{ underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :show-inheritance:

   {% if '__init__' in methods %}
     {% set caught_result = methods.remove('__init__') %}
   {% endif %}

   {% block attributes_summary %}
   {% if attributes %}

   .. rubric:: Attributes Summary

   .. autosummary::
      :toctree:
      :template: base.rst
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}

   {% endif %}
   {% endblock %}

   {% block methods_summary %}
   {% if methods %}

   .. rubric:: Methods Summary

   .. autosummary::
      :toctree:
      :template: base.rst
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}

   {% endif %}
   {% endblock %}
