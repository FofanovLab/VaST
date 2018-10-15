
import pytest
import os
import VaST.analyze


PICK_BEST_SETS_TUPLE = ({'1111':
     {'123': 
          {'primer_zone': {'med_size': 3.5, 'percent_ok': 70.0, 'downstream': "1,0,0,0,1,10,0,0,0,0", 'upstream': "1,0,0,0,1,10,0,0,0,0" },
           's': 123,
           'sites': ['CP020414.2::123::123'],
           'g': {'length': 4280582, 'id': 1, 'name': 'CP020414.2'}},
     '234':
          {'primer_zone': {'med_size': 2.5, 'percent_ok': 70.58823529411765, 'downstream': '1,0,0,0,0,0,0', 'upstream': '0,0,8,4,0,1,10,0,0,0'},
           's': 234,
           'sites': ['CP020414.2::234::234'],
           'g': {'length': 4280582, 'id': 1, 'name': 'CP020414.2'}},
     }},
 {'1111':
     {'123': 
          {'primer_zone': {'med_size': 3.5, 'percent_ok': 70.0, 'downstream': "1,0,0,0,1,10,0,0,0,0", 'upstream': "1,0,0,0,1,10,0,0,0,0" },
           's': 123,
           'sites': ['CP020414.2::123::123'],
           'g': {'length': 4280582, 'id': 1, 'name': 'CP020414.2'}}}},
 {'1111':
     {'123': 
          {'primer_zone': {'med_size': 3.5, 'percent_ok': 70.0, 'downstream': "1,0,0,0,1,10,0,0,0,0", 'upstream': "1,0,0,0,1,10,0,0,0,0" },
           's': 123,
           'sites': ['CP020414.2::123::123'],
           'g': {'length': 4280582, 'id': 1, 'name': 'CP020414.2'}}},
 '2222':
     {'456': 
         {'primer_zone': {'med_size': 4.5, 'percent_ok': 75.0, 'downstream': "1,0,0,0,0,0,0,0,0,0", 'upstream': "1,0,0,0,1,10,0,0,0,0,23,5,1,0,0,0,0,0"},
           's': 456,
           'sites': ['CP020414.2::456::456'],
           'g': {'length': 4280582, 'id': 1, 'name': 'CP020414.2'}}}
     },
 
  {'1111':
     {'123': 
          {'primer_zone': {'med_size': 3.5, 'percent_ok': 70.0, 'downstream': "1,0,0,0,1,10,0,0,0,0", 'upstream': "1,0,0,0,1,10,0,0,0,0" },
           's': 123,
           'sites': ['CP020414.2::123::123'],
           'g': {'length': 4280582, 'id': 1, 'name': 'CP020414.2'}},
     '234':
          {'primer_zone': {'med_size': 2.5, 'percent_ok': 70.58823529411765, 'downstream': '1,0,0,0,0,0,0', 'upstream': '0,0,8,4,0,1,10,0,0,0'},
           's': 234,
           'sites': ['CP020414.2::234::234'],
           'g': {'length': 4280582, 'id': 1, 'name': 'CP020414.2'}}},
 '2222':
     {'456': 
         {'primer_zone': {'med_size': 4.5, 'percent_ok': 75.0, 'downstream': "1,0,0,0,0,0,0,0,0,0", 'upstream': "1,0,0,0,1,10,0,0,0,0,23,5,1,0,0,0,0,0"},
           's': 456,
           'sites': ['CP020414.2::456::456'],
           'g': {'length': 4280582, 'id': 1, 'name': 'CP020414.2'}}}
     })

def test_get_file_from_project_multiple_matches_raises(monkeypatch):
    monkeypatch.setattr(VaST.analyze,
        'expand_files', lambda x, y: ['file1', 'file2'])
    with pytest.raises(ValueError):
        VaST.analyze.get_file_from_project("", "")

def test_get_file_from_project_no_matching_files_raises(monkeypatch):
    monkeypatch.setattr(VaST.analyze,
        'expand_files', lambda x, y: [])
    with pytest.raises(IndexError):
        VaST.analyze.get_file_from_project("", "")

def test_get_file_from_project_single_match_no_file(monkeypatch):
    monkeypatch.setattr(os.path, 'isfile', lambda x: False)
    monkeypatch.setattr(VaST.analyze,
        'expand_files', lambda x, y: ["file"])
    with pytest.raises(FileNotFoundError):
        VaST.analyze.get_file_from_project("", "")

def test_get_file_from_project_single_file_match(monkeypatch):
    monkeypatch.setattr(os.path, 'isfile', lambda x: True)
    monkeypatch.setattr(VaST.analyze,
        'expand_files', lambda x, y: ["file"])
    assert VaST.analyze.get_file_from_project("", "") == "file"

def test_pick_best_loci_empty_dict_return_empty_dict():
    assert VaST.analyze.pick_best_loci({}) == {}

def test_remove_extra_loci_empty_dict_return_empty_dict():
    assert VaST.analyze.remove_extra_loci({}) == {}


@pytest.mark.parametrize("patterns", PICK_BEST_SETS_TUPLE)
def test_remove_extra_loci_return_patterns_with_one_loci(patterns):
    result = VaST.analyze.remove_extra_loci(patterns)
    assert [
        len(loci.keys()) for _, 
        loci in result.items()] == [
            1 for _ in range(len(result))]

