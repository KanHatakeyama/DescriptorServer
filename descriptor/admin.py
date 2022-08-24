from django.contrib import admin
from . import models
# Register your models here.


class MoleculeAdmin(admin.ModelAdmin):
    pass


admin.site.register(models.Molecule, MoleculeAdmin)
