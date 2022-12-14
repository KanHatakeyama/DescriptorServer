# Generated by Django 3.1.2 on 2022-08-24 10:04

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('descriptor', '0003_molecule_avfp'),
    ]

    operations = [
        migrations.RenameField(
            model_name='molecule',
            old_name='avfp',
            new_name='avfp_desc',
        ),
        migrations.AddField(
            model_name='molecule',
            name='jr_desc',
            field=models.TextField(blank=True, max_length=100000, null=True),
        ),
        migrations.AddField(
            model_name='molecule',
            name='mord2d_desc',
            field=models.TextField(blank=True, max_length=100000, null=True),
        ),
    ]
