#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 15:58:32 2026

@author: syedtariqishtiaq
"""

# Learning classes

class Dog:
    def __init__(self, name, age):
        self.name = name
        self.age = age

    def bark(self):
        print(f"{self.name} says woof!")
        
    def sleep(self):
        print(f"{self.name} is tired")
        
rex = Dog(name="Rex", age=3)
bella = Dog(name="Bella", age=5)


bella.bark()
rex.sleep()


class Pet: 
    def __init__(self, name, species):
        self.name = name
        self.species = species
  
dog = Pet(name = 'Fido', species= 'Dog')

'''
Create a class Student:

Attributes: name, grade
Method: display_info() → prints both values
'''

class Student:
    def __init__(self, name, grade):
        self.name = name
        self.grade = grade
        
    def display_info(self):
        print('Name:', self.name,';', 'Grade:', self.grade)
        
s1 = Student('Tariq', 90)
s2 = Student('Ayesha', 80)

s1.display_info()


class Car:
    def __init__(self, brand, speed=0):
        
















