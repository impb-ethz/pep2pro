-- Adminer 3.3.3 MySQL dump

SET NAMES utf8;
SET foreign_key_checks = 0;
SET time_zone = 'SYSTEM';
SET sql_mode = 'NO_AUTO_VALUE_ON_ZERO';

DROP DATABASE IF EXISTS `pep2pro`;
CREATE DATABASE `pep2pro` /*!40100 DEFAULT CHARACTER SET utf8 */;
USE `pep2pro`;

DROP TABLE IF EXISTS `datfile`;
CREATE TABLE `datfile` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `datfile` varchar(255) NOT NULL,
  `ionsscore` double NOT NULL,
  `probability` double NOT NULL,
  `searchdatabase_id` int(11) NOT NULL,
  `experiment` varchar(255) NOT NULL,
  `genotype` varchar(255) NOT NULL,
  `tissue` varchar(255) NOT NULL,
  `treatment` varchar(255) NOT NULL,
  `replicate` varchar(255) NOT NULL,
  `fraction` varchar(255) NOT NULL,
  PRIMARY KEY (`id`),
  KEY `searchdatabase_id` (`searchdatabase_id`),
  CONSTRAINT `datfile_ibfk_1` FOREIGN KEY (`searchdatabase_id`) REFERENCES `searchdatabase` (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;


DROP TABLE IF EXISTS `hit`;
CREATE TABLE `hit` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `datfile_id` int(11) NOT NULL,
  `measurement_id` int(11) NOT NULL,
  `peptide_id` int(11) NOT NULL,
  `peptide_modified_id` int(11) NOT NULL,
  `decoy` bit(1) NOT NULL,
  `ambiguous` bit(1) NOT NULL,
  `multiloci` bit(1) NOT NULL,
  `query` int(11) NOT NULL,
  `rank` int(11) NOT NULL,
  `ionsscore` double NOT NULL,
  `probability` double NOT NULL,
  `scannbr` varchar(100) NOT NULL,
  `charge` int(11) NOT NULL,
  `retention_time` varchar(100) NOT NULL,
  `calcmass` double NOT NULL,
  `deltamass` double NOT NULL,
  `obsmass` double NOT NULL,
  `expmass` double NOT NULL,
  PRIMARY KEY (`id`),
  KEY `datfile_id` (`datfile_id`),
  KEY `measurement_id` (`measurement_id`),
  KEY `peptide_id` (`peptide_id`),
  KEY `peptide_modified_id` (`peptide_modified_id`),
  CONSTRAINT `hit_ibfk_1` FOREIGN KEY (`datfile_id`) REFERENCES `datfile` (`id`) ON DELETE CASCADE,
  CONSTRAINT `hit_ibfk_2` FOREIGN KEY (`measurement_id`) REFERENCES `measurement` (`id`) ON DELETE CASCADE,
  CONSTRAINT `hit_ibfk_3` FOREIGN KEY (`peptide_id`) REFERENCES `peptide` (`id`) ON DELETE CASCADE,
  CONSTRAINT `hit_ibfk_4` FOREIGN KEY (`peptide_modified_id`) REFERENCES `peptide_modified` (`id`) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8;


DROP TABLE IF EXISTS `hit_modificationtype`;
CREATE TABLE `hit_modificationtype` (
  `hit_id` int(11) NOT NULL,
  `modificationtype_id` int(11) NOT NULL,
  KEY `hit_id` (`hit_id`),
  KEY `modificationtype_id` (`modificationtype_id`),
  CONSTRAINT `hit_modificationtype_ibfk_1` FOREIGN KEY (`hit_id`) REFERENCES `hit` (`id`) ON DELETE CASCADE,
  CONSTRAINT `hit_modificationtype_ibfk_2` FOREIGN KEY (`modificationtype_id`) REFERENCES `modificationtype` (`id`) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8;


DROP TABLE IF EXISTS `measurement`;
CREATE TABLE `measurement` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `measurement` varchar(255) NOT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;


DROP TABLE IF EXISTS `modification`;
CREATE TABLE `modification` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `modification` varchar(25) NOT NULL,
  `modificationtype_id` int(11) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `modification` (`modification`),
  KEY `modificationtype_id` (`modificationtype_id`),
  CONSTRAINT `modification_ibfk_1` FOREIGN KEY (`modificationtype_id`) REFERENCES `modificationtype` (`id`) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8;


DROP TABLE IF EXISTS `modificationtype`;
CREATE TABLE `modificationtype` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `modificationtype` varchar(25) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `modificationtype` (`modificationtype`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;


DROP TABLE IF EXISTS `peptide`;
CREATE TABLE `peptide` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `peptide` varchar(255) NOT NULL,
  `missed_cleavages` int(11) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `peptide` (`peptide`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;


DROP TABLE IF EXISTS `peptide_modification`;
CREATE TABLE `peptide_modification` (
  `peptide_modified_id` int(11) NOT NULL,
  `modpos` int(11) NOT NULL,
  `modification_id` int(11) NOT NULL,
  KEY `modification_id` (`modification_id`),
  KEY `peptide_modified_id` (`peptide_modified_id`),
  CONSTRAINT `peptide_modification_ibfk_1` FOREIGN KEY (`modification_id`) REFERENCES `modification` (`id`) ON DELETE CASCADE,
  CONSTRAINT `peptide_modification_ibfk_2` FOREIGN KEY (`peptide_modified_id`) REFERENCES `peptide_modified` (`id`) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8;


DROP TABLE IF EXISTS `peptide_modified`;
CREATE TABLE `peptide_modified` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `peptide_id` int(11) NOT NULL,
  `peptide_modified` varchar(255) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `peptide_modified` (`peptide_modified`),
  KEY `peptide_id` (`peptide_id`),
  CONSTRAINT `peptide_modified_ibfk_2` FOREIGN KEY (`peptide_id`) REFERENCES `peptide` (`id`) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8;


DROP TABLE IF EXISTS `peptide_protein`;
CREATE TABLE `peptide_protein` (
  `peptide_id` int(11) NOT NULL,
  `protein_id` int(11) NOT NULL,
  KEY `peptide_id` (`peptide_id`),
  KEY `protein_id` (`protein_id`),
  CONSTRAINT `peptide_protein_ibfk_1` FOREIGN KEY (`peptide_id`) REFERENCES `peptide` (`id`) ON DELETE CASCADE,
  CONSTRAINT `peptide_protein_ibfk_2` FOREIGN KEY (`protein_id`) REFERENCES `protein` (`id`) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8;


DROP TABLE IF EXISTS `protein`;
CREATE TABLE `protein` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `searchdatabase_id` int(11) NOT NULL,
  `accession` varchar(40) NOT NULL,
  `locus` varchar(40) NOT NULL,
  `description` varchar(255) NOT NULL,
  `tryptic_peptides` int(11) NOT NULL,
  `decoy` bit(1) NOT NULL,
  `protein_length` int(11) NOT NULL,
  `isoelectric_point` double NOT NULL,
  `molecular_weight` double NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `searchdatabase_id_accession` (`searchdatabase_id`,`accession`),
  CONSTRAINT `protein_ibfk_1` FOREIGN KEY (`searchdatabase_id`) REFERENCES `searchdatabase` (`id`) ON DELETE CASCADE
) ENGINE=InnoDB DEFAULT CHARSET=utf8;


DROP TABLE IF EXISTS `searchdatabase`;
CREATE TABLE `searchdatabase` (
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `searchdb` varchar(25) NOT NULL,
  `searchdatabase` varchar(255) NOT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `searchdb` (`searchdb`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

DROP USER 'pep2pro'@'localhost';

CREATE USER 'pep2pro'@'localhost' IDENTIFIED BY 'pep2pro';

GRANT SELECT,INSERT,UPDATE,DELETE,LOCK TABLES ON pep2pro.* TO 'pep2pro'@'localhost';

FLUSH PRIVILEGES;

-- 2018-02-22 14:14:48
