#! /bin/sh
#
# read_flowline_200mhz.sh
# Copyright (C) 2017 davidl <davidl@ryder>
#
# Distributed under terms of the MIT license.
#
gp=0,20,25,25,28,30,30,30,30,30,30,30,30,30

./read_gssi.py --rev 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 --gp $gp --o flowline0-25 --stack 2 ../radar_backup_20170108/SOUTHPOLE__002.DZT ../radar_backup_20170108/SOUTHPOLE__003.DZT ../radar_backup_20170108/SOUTHPOLE__004.DZT ../radar_backup_20170108/SOUTHPOLE__005.DZT

./read_gssi.py --rev 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 --gp $gp --o flowline25-50 --stack 2 ../radar_backup_20170108/SOUTHPOLE__006.DZT ../radar_backup_20170108/SOUTHPOLE__007.DZT ../radar_backup_20170108/SOUTHPOLE__008.DZT

./read_gssi.py --rev 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 --gp $gp --o flowline50-75 --stack 2 ../radar_backup_20170108/SOUTHPOLE__018.DZT ../radar_backup_20170108/SOUTHPOLE__019.DZT

./read_gssi.py --rev 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 --gp $gp --o flowline75-100 --stack 2 ../radar_backup_20170108/SOUTHPOLE__021.DZT ../radar_backup_20170108/SOUTHPOLE__022.DZT

./read_gssi.py --rev 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 --gp $gp --o flowline --stack 10 ../radar_backup_20170108/SOUTHPOLE__002.DZT ../radar_backup_20170108/SOUTHPOLE__003.DZT ../radar_backup_20170108/SOUTHPOLE__004.DZT ../radar_backup_20170108/SOUTHPOLE__005.DZT ../radar_backup_20170108/SOUTHPOLE__006.DZT ../radar_backup_20170108/SOUTHPOLE__007.DZT ../radar_backup_20170108/SOUTHPOLE__008.DZT ../radar_backup_20170108/SOUTHPOLE__018.DZT ../radar_backup_20170108/SOUTHPOLE__019.DZT ../radar_backup_20170108/SOUTHPOLE__021.DZT ../radar_backup_20170108/SOUTHPOLE__022.DZT
