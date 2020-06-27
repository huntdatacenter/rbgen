`edit-bgen` is a tool for editing existing BGEN file metadata.  Currently it supports editing two pieces of data:

* Editing the 'free data area' in the header block.  For example, a use case is stamping BGEN files with a unique code or timestamp.
* Removing sample identifiers from a BGEN file - since sometimes it's preferred to keep sample identifiers seperate.

`edit-bgen` is built along with the rest of the BGEN repo.

**Warning**

`edit-bgen` edits files destructively.  It's therefore possible to lose data using `edit-bgen`; you should backup any file you care about before using it.

