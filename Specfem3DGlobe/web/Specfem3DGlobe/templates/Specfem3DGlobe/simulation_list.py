(
    {% for object in object_list %}
    { 'id': {{ object.id }}, 'status': '{{ object.status }}' },
    {% endfor %}
)
