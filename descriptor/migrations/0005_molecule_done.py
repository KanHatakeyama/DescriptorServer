# Generated by Django 3.1.2 on 2022-08-26 04:07

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('descriptor', '0004_auto_20220824_1004'),
    ]

    operations = [
        migrations.AddField(
            model_name='molecule',
            name='done',
            field=models.BooleanField(default=False),
        ),
    ]
