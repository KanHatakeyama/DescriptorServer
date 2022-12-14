from django.db import models

# Create your models here.


class Molecule(models.Model):
    SMILES = models.CharField(max_length=1500)

    RDKit_desc = models.TextField(null=True, blank=True, max_length=10**5)
    avfp_desc = models.TextField(null=True, blank=True, max_length=10**5)
    mord2d_desc = models.TextField(null=True, blank=True, max_length=10**5)
    jr_desc = models.TextField(null=True, blank=True, max_length=10**5)
    PM7_desc = models.TextField(
        null=True, blank=True, max_length=10**5)

    done = models.BooleanField(default=False)
    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)

    def __str__(self) -> str:
        return self.SMILES
