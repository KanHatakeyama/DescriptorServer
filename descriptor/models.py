from django.db import models

# Create your models here.


class Molecule(models.Model):
    SMILES = models.CharField(max_length=1500)
    RDKit_desc = models.TextField(null=True, blank=True, max_length=10**5)

    def __str__(self) -> str:
        return self.SMILES
