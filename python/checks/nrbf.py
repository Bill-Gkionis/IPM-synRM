"""nrbf: read a .NET BinaryFormatter file (MS-NRBF) without starting any .NET program. Used to read PrePoMax's settings.bin.

What it does, and why
    PrePoMax keeps its settings (which solver it runs, which work folder it uses) in a file called settings.bin. The file is the
    output of .NET's BinaryFormatter: a stream of records that describe objects, their classes, their member names and their
    values. The brief wants to know which solver PrePoMax would run (the `_calculixExe` setting) WITHOUT starting PrePoMax, whose
    window would open on the desktop. So this module reads the stream itself. The format is Microsoft's published
    specification [MS-NRBF]; the reading below follows it record by record.

The format, as far as it is needed (every number is little-endian)
    A stream is a sequence of records. Each starts with one byte, the record type:
        0  SerializedStreamHeader   root id, header id, major version, minor version (four Int32): must be the first record
        12 BinaryLibrary            library id (Int32), the library's name (a string): says where the classes come from
        5  ClassWithMembersAndTypes object id, class name, member count, member names; then a type for each member; then the
                                    library id; then the values. (4 is the same without a library id, for system classes.)
        1  ClassWithId              object id and the id of an earlier class: the same members, new values (the cheap repeat)
        6  BinaryObjectString       object id, then the string
        9  MemberReference          the id of an object defined somewhere in the stream (a repeated string is written once)
        10 ObjectNull, 13/14 ObjectNullMultiple   one or several nulls
        8  MemberPrimitiveTyped     a type code and a value (a boxed number)
        7  BinaryArray, 15 ArraySinglePrimitive, 16 ArraySingleObject, 17 ArraySingleString
        11 MessageEnd               the last record
    A string is written as its length in bytes, in 7 bits per byte with the top bit meaning "more follows", then its UTF-8 bytes.
    A member's value is inline for a primitive type (an Int32 is four bytes) and a record for everything else.

What it checks (a reader that guesses is worse than none)
    It reads to MessageEnd, and then the stream must be exactly used up: not a byte missing, not a byte left. Any record it does
    not know, any length that runs past the end, any reference to an object that was never defined is an NrbfError that says what
    and where (the byte offset). A caller therefore gets either the whole file read correctly or an error with a reason.

Self-check (the answers follow).
  1. The string "ccx" is written as 03 63 63 78. What is the string of 200 bytes written as in its first bytes?
  2. Why is a MemberReference needed?
  3. Why does the reader insist that no byte is left after MessageEnd?
Answers.
  1. C8 01 followed by the 200 bytes: 200 = 0b1_1001000 is 7 bits (1001000 = 0x48) with the top bit set to say "more follows"
     (0xC8) and then the next 7 bits (0000001 = 0x01).
  2. The same string object can be a member of many objects; the formatter writes it once and refers to it by id afterwards.
     A reader that only handled inline strings would see a number where it expected text.
  3. A leftover byte means the reader lost step with the writer somewhere (a record it understood only partly). The first
     member values may still look right, so the only safe test of a whole-file read is that it ends exactly where the file does.
"""
from __future__ import annotations

import struct

__all__ = ["NrbfError", "read_class_members"]


class NrbfError(Exception):
    """The stream is not a BinaryFormatter stream that this reader can read completely; the message says why and where."""


class Obj:
    """A class instance of the stream: its class name and its members, in order."""

    def __init__(self, ident: int, cls: str, names: list[str]):
        self.ident, self.cls, self.names, self.values = ident, cls, names, []

    def __repr__(self) -> str:
        return f"Obj({self.cls}, {dict(zip(self.names, self.values))})"


class Ref:
    def __init__(self, ident: int):
        self.ident = ident


class Null:
    def __init__(self, count: int = 1):
        self.count = count


_FORMATS = {1: ("?", 1), 2: ("B", 1), 6: ("<d", 8), 7: ("<h", 2), 8: ("<i", 4), 9: ("<q", 8), 10: ("b", 1), 11: ("<f", 4),
            12: ("<q", 8), 13: ("<Q", 8), 14: ("<H", 2), 15: ("<I", 4), 16: ("<Q", 8)}          # primitive type code -> (struct format, size)


class _Reader:
    def __init__(self, data: bytes):
        self.data, self.pos = data, 0
        self.objects: dict[int, object] = {}            # object id -> Obj | str | list
        self.classes: dict[int, tuple[str, list[str], list]] = {}      # class id -> (name, member names, member types)
        self.instances: list[Obj] = []

    def fail(self, why: str) -> NrbfError:
        return NrbfError(f"{why} (byte {self.pos} of {len(self.data)})")

    def take(self, n: int) -> bytes:
        if n < 0 or self.pos + n > len(self.data):
            raise self.fail(f"the stream ends in the middle of a record: {n} more byte(s) needed")
        chunk = self.data[self.pos:self.pos + n]
        self.pos += n
        return chunk

    def u8(self) -> int:
        return self.take(1)[0]

    def i32(self) -> int:
        return struct.unpack("<i", self.take(4))[0]

    def string(self) -> str:
        length, shift = 0, 0
        while True:
            byte = self.u8()
            length |= (byte & 0x7F) << shift
            if not byte & 0x80:
                break
            shift += 7
            if shift > 28:
                raise self.fail("a string length that is not a valid 7-bit number")
        try:
            return self.take(length).decode("utf-8")
        except UnicodeDecodeError:
            raise self.fail("a string that is not UTF-8") from None

    def primitive(self, code: int):
        if code in (5, 18):                              # Decimal and String are written as strings
            return self.string()
        if code == 3:                                    # Char: one UTF-8 encoded character, 1 to 4 bytes
            first = self.u8()
            extra = 0 if first < 0x80 else 1 if first < 0xE0 else 2 if first < 0xF0 else 3
            self.pos -= 1
            return self.take(1 + extra).decode("utf-8")
        if code not in _FORMATS:
            raise self.fail(f"unknown primitive type {code}")
        fmt, size = _FORMATS[code]
        return struct.unpack(fmt, self.take(size))[0]

    def class_info(self) -> tuple[int, str, list[str]]:
        ident, name, count = self.i32(), self.string(), self.i32()
        if not 0 <= count <= 100000:
            raise self.fail(f"{count} members in class {name!r}")
        return ident, name, [self.string() for _ in range(count)]

    def additional(self, kind: int):
        """The extra information a binary type carries: a primitive code, a class name, or a class name and library."""
        if kind in (0, 7):
            return self.u8()
        if kind == 3:
            return self.string()
        if kind == 4:
            return (self.string(), self.i32())
        return None

    def member_types(self, count: int) -> list[tuple[int, object]]:
        kinds = [self.u8() for _ in range(count)]
        return [(kind, self.additional(kind)) for kind in kinds]

    def members(self, obj: Obj, types: list[tuple[int, object]]) -> Obj:
        self.objects[obj.ident] = obj
        self.instances.append(obj)
        pending_nulls = 0
        for kind, extra in types:
            if pending_nulls:                            # the rest of an ObjectNullMultiple that covers several members
                pending_nulls -= 1
                obj.values.append(None)
                continue
            if kind == 0:
                obj.values.append(self.primitive(extra))
                continue
            value = self.record()
            if isinstance(value, Null):
                pending_nulls = value.count - 1
                value = None
            obj.values.append(value)
        return obj

    def record(self):
        kind = self.u8()
        if kind == 0:
            self.take(16)                                # root id, header id, major and minor version
            return self.record()
        if kind == 12:                                   # BinaryLibrary: id and name, then the record it is a prefix of
            self.i32()
            self.string()
            return self.record()
        if kind in (4, 5):
            ident, name, names = self.class_info()
            types = self.member_types(len(names))
            if kind == 5:
                self.i32()                               # the library id
            self.classes[ident] = (name, names, types)
            return self.members(Obj(ident, name, names), types)
        if kind == 1:
            ident, metadata = self.i32(), self.i32()
            if metadata not in self.classes:
                raise self.fail(f"ClassWithId refers to class {metadata}, which was not defined")
            name, names, types = self.classes[metadata]
            return self.members(Obj(ident, name, names), types)
        if kind == 6:
            ident = self.i32()
            text = self.string()
            self.objects[ident] = text
            return text
        if kind == 8:
            return self.primitive(self.u8())
        if kind == 9:
            return Ref(self.i32())
        if kind == 10:
            return Null()
        if kind == 13:
            return Null(self.u8())
        if kind == 14:
            return Null(self.i32())
        if kind == 15:                                   # ArraySinglePrimitive: values inline
            ident, length, code = self.i32(), self.i32(), self.u8()
            items = [self.primitive(code) for _ in range(length)]
            self.objects[ident] = items
            return items
        if kind in (16, 17):                             # ArraySingleObject / ArraySingleString: records
            ident, length = self.i32(), self.i32()
            items = self.array_items(length)
            self.objects[ident] = items
            return items
        if kind == 7:                                    # BinaryArray
            ident, array_type, rank = self.i32(), self.u8(), self.i32()
            lengths = [self.i32() for _ in range(rank)]
            if array_type in (3, 4, 5):                  # the array types with lower bounds
                for _ in range(rank):
                    self.i32()
            element_kind = self.u8()
            extra = self.additional(element_kind)
            total = 1
            for length in lengths:
                total *= length
            if element_kind == 0:
                items = [self.primitive(extra) for _ in range(total)]
            else:
                items = self.array_items(total)
            self.objects[ident] = items
            return items
        raise self.fail(f"unknown record type {kind}")

    def array_items(self, length: int) -> list:
        items: list = []
        while len(items) < length:
            value = self.record()
            if isinstance(value, Null):
                if len(items) + value.count > length:        # checked before the list is made: a count of 2 billion must not allocate
                    raise self.fail("a run of nulls goes past the length of its array")
                items.extend([None] * value.count)
            else:
                items.append(value)
        return items

    def read_all(self) -> None:
        if not self.data or self.data[0] != 0:
            raise NrbfError("the stream does not start with a SerializedStreamHeader (is it a BinaryFormatter file?)")
        while True:
            if self.pos >= len(self.data):
                raise self.fail("the stream ends before MessageEnd")
            if self.data[self.pos] == 11:
                self.pos += 1
                break
            self.record()
        if self.pos != len(self.data):
            raise self.fail(f"{len(self.data) - self.pos} byte(s) left after MessageEnd")

    def resolve(self, value):
        """A member value as plain Python: a reference becomes the object it names (a string, an Obj or a list)."""
        if isinstance(value, Ref):
            if value.ident not in self.objects:
                raise NrbfError(f"a reference to object {value.ident}, which the stream never defines")
            return self.objects[value.ident]
        return value


def read_class_members(data: bytes, class_name: str) -> dict[str, object]:
    """{member name: value} of the first object of class `class_name` in a BinaryFormatter stream.

    A primitive member is a Python bool, int or float; a string member (inline or by reference) is a str; a null is None; an
    object member is the Obj itself. NrbfError, with a reason, if the stream is not read completely or has no such class."""
    reader = _Reader(data)
    try:
        reader.read_all()
    except NrbfError:
        raise
    except (struct.error, ValueError, KeyError, IndexError, RecursionError) as error:
        raise reader.fail(f"cannot read the stream: {type(error).__name__}: {error}") from None
    for obj in reader.instances:
        if obj.cls == class_name:
            return {name: reader.resolve(value) for name, value in zip(obj.names, obj.values)}
    raise NrbfError(f"the stream has no object of class {class_name!r} (it has: {', '.join(sorted({o.cls for o in reader.instances})) or 'none'})")
