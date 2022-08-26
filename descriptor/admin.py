from django.contrib import admin
from . import models
# Register your models here.


class MoleculeAdmin(admin.ModelAdmin):
    list_filter = ["done",
                   "created", "updated"]

    list_display = ["SMILES", "done",
                    "created", "updated"]


admin.site.register(models.Molecule, MoleculeAdmin)
