// Minimal, dependency-free ZIP archive builder (store method, no compression).
//
// Figures are already compact (PNG is compressed; a handful of SVGs are small),
// so storing them uncompressed keeps the writer simple and bulletproof. Used to
// bundle "Save all figures" into one download instead of N blocked downloads.

const CRC_TABLE = (() => {
    const table = new Uint32Array(256);
    for (let n = 0; n < 256; n += 1) {
        let c = n;
        for (let k = 0; k < 8; k += 1) {
            c = (c & 1) ? (0xedb88320 ^ (c >>> 1)) : (c >>> 1);
        }
        table[n] = c >>> 0;
    }
    return table;
})();

const crc32 = (bytes) => {
    let crc = 0xffffffff;
    for (let i = 0; i < bytes.length; i += 1) {
        crc = CRC_TABLE[(crc ^ bytes[i]) & 0xff] ^ (crc >>> 8);
    }
    return (crc ^ 0xffffffff) >>> 0;
};

const encoder = new TextEncoder();

// entries: [{ name: string, data: Uint8Array }] -> a Blob of a .zip file.
export const buildZip = (entries) => {
    const files = entries.map((entry) => ({
        nameBytes: encoder.encode(entry.name),
        data: entry.data,
        crc: crc32(entry.data),
    }));

    const LOCAL_HEADER = 30;
    const CENTRAL_HEADER = 46;
    const EOCD = 22;
    let size = EOCD;
    files.forEach((file) => {
        size += LOCAL_HEADER + file.nameBytes.length + file.data.length;
        size += CENTRAL_HEADER + file.nameBytes.length;
    });

    const buffer = new ArrayBuffer(size);
    const view = new DataView(buffer);
    const out = new Uint8Array(buffer);
    let offset = 0;
    const localOffsets = [];

    const u32 = (value) => { view.setUint32(offset, value, true); offset += 4; };
    const u16 = (value) => { view.setUint16(offset, value, true); offset += 2; };
    const bytes = (data) => { out.set(data, offset); offset += data.length; };

    files.forEach((file) => {
        localOffsets.push(offset);
        u32(0x04034b50); // local file header signature
        u16(20);         // version needed to extract
        u16(0);          // general purpose flags
        u16(0);          // compression method: store
        u16(0);          // last mod time
        u16(0);          // last mod date
        u32(file.crc);
        u32(file.data.length); // compressed size
        u32(file.data.length); // uncompressed size
        u16(file.nameBytes.length);
        u16(0);          // extra field length
        bytes(file.nameBytes);
        bytes(file.data);
    });

    const centralStart = offset;
    files.forEach((file, index) => {
        u32(0x02014b50); // central directory header signature
        u16(20);         // version made by
        u16(20);         // version needed to extract
        u16(0);          // general purpose flags
        u16(0);          // compression method
        u16(0);          // last mod time
        u16(0);          // last mod date
        u32(file.crc);
        u32(file.data.length);
        u32(file.data.length);
        u16(file.nameBytes.length);
        u16(0);          // extra field length
        u16(0);          // file comment length
        u16(0);          // disk number start
        u16(0);          // internal file attributes
        u32(0);          // external file attributes
        u32(localOffsets[index]);
        bytes(file.nameBytes);
    });
    const centralSize = offset - centralStart;

    u32(0x06054b50); // end of central directory signature
    u16(0);          // number of this disk
    u16(0);          // disk with central directory
    u16(files.length);
    u16(files.length);
    u32(centralSize);
    u32(centralStart);
    u16(0);          // comment length

    return new Blob([buffer], { type: 'application/zip' });
};
