(
    {% for object in object_list %}
    { 'id': {{ object.id }}, 'status': '{{ object.get_status_display }}' },
    {% endfor %}
)
