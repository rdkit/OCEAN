from django.contrib import admin
from django.db import models
from ocean.models import *


class DataSourcesAdmin(admin.ModelAdmin):

    list_display = ('id','name')

admin.site.register(DataSources,DataSourcesAdmin)

if os.path.exists('ocean/custom_admin.py'):
    import custom_admin
    mcu = {k:v for k,v in custom_admin.__dict__.items() if not k.startswith('__')}
    current_locals = locals()
    for entry,value in mcu.items():
        if not entry in current_locals.keys() or \
                        entry in []:
            current_locals.update({entry:value})
            print "load custom admin",entry,value
