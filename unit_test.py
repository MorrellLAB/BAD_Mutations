#!/usr/bin/env python
import genome_jgi
import unittest
import os

class TestJGIFunctions(unittest.TestCase):

    def setUp(self):
        self.jgi_user = os.environ['jgi_user']
        self.jgi_password = os.environ['jgi_password']
        self.jgi = genome_jgi.JGIUtils(self.jgi_user,self.jgi_password)

    def test_user(self):
        self.assertTrue( self.jgi_user != None)
        self.assertTrue( self.jgi_password != None)

    def test_fetch_xml(self):
        xml = self.jgi.fetch_xml()
        self.assertTrue( xml != None )
        self.assertTrue( xml.find("cds") )

    def test_url_list(self):
        url_list = self.jgi.fetch_url_list()
        self.assertTrue( len(url_list) > 0 )

    def test_cds_list(self):
        cds_list = self.jgi.fetch_cds_list()
        self.assertTrue( len(cds_list) > 0 )
        for url in cds_list:
            part = url.split(".")
            self.assertTrue('cds' in part)

if __name__ == '__main__':
      unittest.main()
